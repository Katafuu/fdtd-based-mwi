% tds_2d_circular Circular-array transmission-based microwave imaging demo.
%
% This script keeps the same FDTD setup as the earlier planar script, but
% replaces the separate TX/RX line arrays with one circular antenna array.
% Each antenna acts as a transmitter in turn; all other antennas act as
% receivers for that run.
%
% Internal matrix convention used by this project:
%   Ez(x,y,t)
%   mask4D(x,y,tx,rx)
%   epsr_final(x,y,tx)
% For plotting, use plotXY(...), which transposes for MATLAB display only.

%% Repository / path setup
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end

algorithmDir = scriptDir;
workspaceRoot = fileparts(fileparts(fileparts(algorithmDir)));

addpath(fullfile(workspaceRoot, 'buildLib'));
addpath(fullfile(workspaceRoot, 'forward_solver', 'mex'));
addpath(fullfile(algorithmDir, 'lib'));

resultsDir = fullfile(algorithmDir, 'results'); %#ok<NASGU>
inputDir = fullfile(resultsDir, 'mats'); %#ok<NASGU>

assert(isfile(fullfile(workspaceRoot, 'forward_solver', 'Makefile')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'fdtd')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'utility')), ...
       'Could not find forward_solver at the workspace root.');

eps0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c0 = 1/sqrt(mu0*eps0);

%% Simulation parameters
Nx = 400;
Ny = 400;
dx = 0.005;
dy = dx;
Nt = 600;
sourceFreq = 1.5e9;
dt = 0.99/(c0*sqrt(1/dx^2 + 1/dy^2));

ax = 1.0;
ay = 1.0;
az = 1.0;

background = struct("epsr", 1.0, "murx", 1.0, "mury", 1.0, ...
    "cond_e", 0.0, "cond_m", 0.0);

grid = fdtdmat.createGrid([Nx Ny], [dx dy], background, [0 0]);
highGrid = fdtdmat.createGrid(2 .* [Nx Ny], [dx/2 dy/2], background, [0 0]);

% Test target: centered circular inclusion
targetCenter = [(Nx - 1+40) * dx / 2, (Ny - +30) * dy / 2];
targetRadius = 30 * dx;
target = fdtdgeom.shape_circle(highGrid, targetCenter, targetRadius, "physical");
highGrid = fdtdmat.applyRegion(highGrid, target, struct("epsr", 5.0, "cond_e", 0.0));

ezRows = 1:2:(2*Nx - 1);
ezCols = 1:2:(2*Ny - 1);
hxRows = 1:2:(2*Nx - 1);
hxCols = 2:2:(2*Ny - 2);
hyRows = 2:2:(2*Nx - 2);
hyCols = 1:2:(2*Ny - 1);

grid.epsr = highGrid.epsr(ezRows, ezCols);
grid.cond_e = highGrid.cond_e(ezRows, ezCols);
grid.cond_m = highGrid.cond_m(ezRows, ezCols);
grid.murx = highGrid.murx(hxRows, hxCols);
grid.mury = highGrid.mury(hyRows, hyCols);

%% PML construction
pmlThicknessRatio = 0.25;
pmlThickness = floor(Ny * pmlThicknessRatio);
pmlOrder = 6;
targetReflection = 1e-60;
kappaMax = 6.0;
etaBackground = sqrt((mu0 * grid.background.murx) / (eps0 * grid.background.epsr));
pmlPhysicalThickness = pmlThickness*dx;
sigmaMax = fdtdpml.utility_sigmaMaxFromReflection(pmlOrder, ...
    targetReflection, etaBackground, pmlPhysicalThickness);
pml = fdtdpml.build_rectangularPML(grid, pmlThickness, pmlOrder, sigmaMax, kappaMax);

%% Config struct for MEX solver
cfg = struct();
cfg.c0 = c0;
cfg.Nx = Nx;
cfg.Ny = Ny;
cfg.sizeZ = 1;
cfg.Nt = Nt;
cfg.dx = dx;
cfg.dy = dy;
cfg.dt = dt;
cfg.snapshotStart = 0;
cfg.snapshotStride = 1;
cfg.returnEz = true;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = false;

cfg.grid = grid;

cfg.pml = pml;
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = ax;
cfg.pml.ay = ay;
cfg.pml.az = az;

disp("Built FDTD case in workspace: cfg, grid, pml, Nx, Ny, dx, dy, dt.");

%% Circular antenna array
% One shared antenna list. Antenna tx acts as transmitter; all other
% antennas are receivers for that tx.
Nant = 12;
antennaInsetCells = 2;

antenna_center = [round((Nx + 1)/2), round((Ny + 1)/2)];
xMin = pmlThickness + 1 + antennaInsetCells;
xMax = Nx - pmlThickness - antennaInsetCells;
yMin = pmlThickness + 1 + antennaInsetCells;
yMax = Ny - pmlThickness - antennaInsetCells;
antenna_radius_cells = floor(min([antenna_center(1)-xMin, ...
    xMax-antenna_center(1), antenna_center(2)-yMin, ...
    yMax-antenna_center(2)]));
cfg.antennas = struct('numAntennas', Nant, ...
    'center', antenna_center, 'radius', antenna_radius_cells, 'focusPadding', 0);
cfg = buildCircularAntennaArrayIdx(cfg);
antenna_cell_indices = cfg.antennas.pos;

fprintf('Built circular array: Nant = %d, radius = %.1f cells.\n', ...
    Nant, antenna_radius_cells);

figure;
plotAntennaArray(antenna_cell_indices, antenna_center, Nx, Ny, pmlThickness);
title('Circular antenna array');

% Sample the source once; each scan activates only the current transmitter row.
cfg.source = struct();
cfg.source.frequency = sourceFreq;
cfg.source.location = cfg.antennas.pos;
cfg.source.func = @(physicalTime) ...
    (1 - 2 .* (pi .* (sourceFreq .* physicalTime - 1)).^2) .* ...
    exp(-(pi .* (sourceFreq .* physicalTime - 1)).^2);
sourceWaveform = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);
cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);

%% Pair geometry: distance, global ray angle, and relative off-boresight angles
% distance(tx,rx): physical distance from antenna tx to antenna rx
% beta_tx(tx,rx): angle between tx boresight and outgoing ray tx->rx
% beta_rx(tx,rx): angle between rx boresight and incoming ray rx->tx
% pairWeight(tx,rx): cos(beta_tx)*cos(beta_rx), with nonpositive values removed
[distance, globalAngle, beta_tx, beta_rx, pairWeight] = buildCircularPairGeometry( ...
    antenna_cell_indices, dx, dy, antenna_center);

maxAngle = deg2rad(40);
pairAccepted = abs(beta_tx) <= maxAngle & abs(beta_rx) <= maxAngle & isfinite(distance);
pairAccepted(eye(Nant) == 1) = false;

% Optional diagnostics
fprintf('Accepted TX-RX pairs: %d out of %d possible off-diagonal pairs.\n', ...
    nnz(pairAccepted), Nant*(Nant-1));

%% Run object/reference simulations and extract excess delay steps
delaysteps = nan(Nant, Nant);
arrival_obj_steps = nan(Nant, Nant);
arrival_inc_steps = nan(Nant, Nant);

for tx = 1:Nant
    cfg.source.samples(:) = 0;
    cfg.source.samples(tx, :) = sourceWaveform;

    rxList = setdiff(1:Nant, tx);
    receiver_indices_this_tx = antenna_cell_indices(rxList, :);

    % Object scan
    fdtdMexResult = fdtd_mex(cfg);

    % Incident/reference scan in background medium
    incident_cfg = fdtdmat.setBackgroundDefault(cfg, background);
    fdtdMexResult_inc = fdtd_mex(incident_cfg);

    % Parse Ez into internal convention Ez(x,y,t)
    Ez = fdtdMexResult.Ez;

    Ez_inc = fdtdMexResult_inc.Ez;

    arrival_obj = getArrivalSteps(Ez, receiver_indices_this_tx, cfg);
    arrival_inc = getArrivalSteps(Ez_inc, receiver_indices_this_tx, incident_cfg);

    arrival_obj_steps(tx, rxList) = arrival_obj.';
    arrival_inc_steps(tx, rxList) = arrival_inc.';
    delaysteps(tx, rxList) = arrival_obj.' - arrival_inc.';

    fprintf('Finished transmitter %d / %d.\n', tx, Nant);
end

%% Compute pairwise average permittivity estimates
% delaysteps is already arrival_obj - arrival_inc in units of time steps.
% Therefore do NOT subtract 1 here.
Delta_t = delaysteps * dt;

clampNegativeDelays = true;
if clampNegativeDelays
    Delta_t(Delta_t < 0) = 0;
end

eps_r = (1 + (c0 .* Delta_t) ./ distance).^2;
eps_r(eye(Nant) == 1) = NaN;

%% Build paper-style footprint masks for the circular array
xRange = pmlThickness+1 : Nx-pmlThickness;
yRange = pmlThickness+1 : Ny-pmlThickness;

% Use footprint size based on circular antenna arc spacing as a reasonable
% first value. Tune Lfp_cells depending on desired coverage/resolution.
arcSpacingCells = 2*pi*antenna_radius_cells / Nant;
Lfp_cells = max(1, round(arcSpacingCells/2));

[mask4D, xRange, yRange] = buildFootprintMask( ...
    Nx, Ny, pmlThickness, antenna_cell_indices, antenna_cell_indices, Lfp_cells);

% Remove self-pair masks.
for ant = 1:Nant
    mask4D(:,:,ant,ant) = false;
end

fprintf('Built circular footprint masks with Lfp_cells = %d.\n', Lfp_cells);

% Optional midpoint/coverage diagnostic
figure;
plotCircularMidpoints(antenna_cell_indices, pairAccepted, Nx, Ny, pmlThickness);
title('Accepted TX-RX midpoint/footprint centers');

%% Reconstruct weighted epsr maps
[epsr_final, epsr_num, epsr_den] = reconstructEpsrFromWeights( ...
    eps_r, pairWeight, pairAccepted, mask4D);

epsr_avg = averageEpsrFinal(epsr_num, epsr_den);

%% Pad non-PML reconstruction back to full grid size for visualization
epsr_avg_full = ones(Nx, Ny);
epsr_avg_full(xRange, yRange) = epsr_avg;

figure;
plotXY(epsr_avg_full, 1:Nx, 1:Ny);
title('Circular-array weighted average \epsilon_r reconstruction, full grid');

%% Optional verification: constant pair estimates should reconstruct constant values where covered
% eps_r_test = ones(size(eps_r));
% [~, test_num, test_den] = reconstructEpsrFromWeights(eps_r_test, pairWeight, pairAccepted, mask4D);
% test_avg = averageEpsrFinal(test_num, test_den);
% disp([min(test_avg(:), [], 'omitnan'), max(test_avg(:), [], 'omitnan')]);

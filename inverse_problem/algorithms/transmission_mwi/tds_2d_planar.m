% Auto-generated plain MATLAB script from tds_2d.mlx
% Contains the modified TDS reconstruction pipeline.

% build_fdtd_case Build the shared TMz FDTD case in memory.
%
% This script creates the scalar simulation parameters, material maps, PML
% maps, and a MEX-ready cfg struct. It intentionally performs no file I/O and
% does not run either solver front end.

%% ------------------------------------------------------------------------

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

resultsDir = fullfile(algorithmDir, 'results');
inputDir = fullfile(resultsDir, 'mats');

assert(isfile(fullfile(workspaceRoot, 'forward_solver', 'Makefile')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'fdtd')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'utility')), ...
       'Could not find forward_solver at the workspace root.');

eps0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c0 = 1/sqrt(mu0*eps0);

%% ------------------------------------------------------------------------

Nx = 400;
Ny = 400;
dx = 0.005; % sweet spot for freq of 1.5GHz to be represented by 20 points per wavelength, aka physical domain bigger, less fine sampling of wave
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

targetCenter = [(Nx - 1) * dx / 2, (Ny - 1) * dy / 2];
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

%% ------------------------------------------------------------------------

pmlThicknessRatio = 0.25;

pmlThickness = floor(Ny * pmlThicknessRatio);
pmlOrder = 6;
targetReflection = 1e-60;
kappaMax = 6.0;
etaBackground = sqrt((mu0 * grid.background.murx) / (eps0 * grid.background.epsr));
pmlCenter = [(Nx + 1)/2 (Ny + 1)/2];

pmlPhysicalThickness = pmlThickness*dx;
sigmaMax = fdtdpml.utility_sigmaMaxFromReflection(pmlOrder, ...
    targetReflection, etaBackground, pmlPhysicalThickness);
% pml = fdtdpml.build_radialPML(grid, pmlCenter, pmlThickness, ...
    % pmlOrder, sigmaMax, kappaMax, "index");
pml = fdtdpml.build_rectangularPML(grid, pmlThickness, pmlOrder, sigmaMax, kappaMax);

%% ------------------------------------------------------------------------

cfg = struct();
cfg.c0 = c0;
cfg.Nx = Nx;
cfg.Ny = Ny;
cfg.sizeZ = 1;
cfg.Nt = Nt;
cfg.dx = dx;
cfg.dy = dy;
cfg.dt = dt;
cfg.snapshotStart = 0; % which 
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

%% ------------------------------------------------------------------------

N = 50;
Ntx = N;
Nrx = N;
paddingx = 0;
paddingy = 0;

xMin = pmlThickness + 1 + paddingx;
xMax = Nx - pmlThickness - paddingx;

yMin = pmlThickness + 1 + paddingy;
yMax = Ny - pmlThickness - paddingy;

step_size = floor((yMax - yMin + 1) / N);

transmitter_cell_indices = zeros(N, 2);
receiver_cell_indices    = zeros(N, 2);

% Antennas placed evenly inside the non-PML region
y_positions = round(yMin + step_size/2 + (0:N-1)*step_size);

% Left side transmitters
transmitter_cell_indices(:,1) = xMin;
transmitter_cell_indices(:,2) = y_positions(:);

% Right side receivers
receiver_cell_indices(:,1) = xMax;
receiver_cell_indices(:,2) = y_positions(:);

% Sample the source once; each scan activates only the current transmitter row.
cfg.antennas = struct();
cfg.antennas.numAntennas = Ntx;
cfg.antennas.pos = transmitter_cell_indices;
cfg.source = struct();
cfg.source.frequency = sourceFreq;
cfg.source.location = cfg.antennas.pos;
cfg.source.func = @(physicalTime) ...
    (1 - 2 .* (pi .* (sourceFreq .* physicalTime - 1)).^2) .* ...
    exp(-(pi .* (sourceFreq .* physicalTime - 1)).^2);
sourceWaveform = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);
cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);

 % Building distance matrix
% Convert index positions to physical coordinates
tx_x = (transmitter_cell_indices(:,1) - 1) * dx;
tx_y = (transmitter_cell_indices(:,2) - 1) * dy;

rx_x = (receiver_cell_indices(:,1) - 1) * dx;
rx_y = (receiver_cell_indices(:,2) - 1) * dy;

% distance(i,j) = distance from transmitter i to receiver j
distance = sqrt((tx_x - rx_x.').^2 + (tx_y - rx_y.').^2);

delta_x = rx_x.' - tx_x;
delta_y = rx_y.' - tx_y;

angle = atan2(delta_y, delta_x);

%% ------------------------------------------------------------------------

delaysteps = nan(Ntx, Nrx);
tic;
for tx = 1:Ntx
    cfg.source.samples(:) = 0;
    cfg.source.samples(tx, :) = sourceWaveform;

    fdtdMexResult = fdtd_mex(cfg);

    incident_cfg = fdtdmat.setBackgroundDefault(cfg, background);
    fdtdMexResult_inc = fdtd_mex(incident_cfg);

    % parse Ez as you already do
    Ez = fdtdMexResult.Ez;
    Ez_inc = fdtdMexResult_inc.Ez;
    arrival_obj = getArrivalSteps(Ez, receiver_cell_indices, cfg);
    arrival_inc = getArrivalSteps(Ez_inc, receiver_cell_indices, incident_cfg);
    delaysteps(tx,:) = arrival_obj - arrival_inc;
end
simRunTime = toc;
%% ------------------------------------------------------------------------
tic;
% delaysteps is already arrival_obj - arrival_inc in units of time steps.
% Therefore do NOT subtract 1 here.
Delta_t = delaysteps * dt;

% Optional numerical cleanup: for an air/reference scan, small negative
% values can appear from discretization. Leave this off if you want the raw
% measured delay exactly as-is.
clampNegativeDelays = true;
if clampNegativeDelays
    Delta_t(Delta_t < 0) = 0;
end

eps_r = (1 + (c0 .* Delta_t) ./ distance).^2;

%% ------------------------------------------------------------------------

% ------------------------------------------------------------
% Construct mask4D(x, y, transmitter_idx, receiver_idx)
% ------------------------------------------------------------
%
% The TDS papers map each TX-RX estimate to a localized footprint centered
% at the midpoint/intersection point of that TX-RX path. A full thick-ray
% mask is useful for debugging, but it produces fan/cone smearing artifacts.

Ntx = size(transmitter_cell_indices, 1);
Nrx = size(receiver_cell_indices, 1);

usePaperFootprintMask = true;

% Footprint half-width in cells. Tune this parameter.
% The paper uses physical footprint sizes such as 7 mm or 9 mm; here we use
% a grid-cell equivalent based on antenna spacing as a first practical value.
Lfp_cells = max(1, round(step_size/2));

if usePaperFootprintMask
    [mask4D, xRange, yRange] = buildFootprintMask( ...
        Nx, Ny, pmlThickness, ...
        transmitter_cell_indices, receiver_cell_indices, Lfp_cells);
else
    % Optional full thick-ray mask for comparison/debugging only.
    xRange = pmlThickness+1 : Nx-pmlThickness;
    yRange = pmlThickness+1 : Ny-pmlThickness;
    [Xidx, Yidx] = ndgrid(xRange, yRange);

    rayWidth_cells = max(step_size, 1);
    halfWidth = rayWidth_cells / 2;

    mask4D = false(numel(xRange), numel(yRange), Ntx, Nrx);

    for tx = 1:Ntx
        x1 = transmitter_cell_indices(tx,1);
        y1 = transmitter_cell_indices(tx,2);

        for rx = 1:Nrx
            x2 = receiver_cell_indices(rx,1);
            y2 = receiver_cell_indices(rx,2);

            vx = x2 - x1;
            vy = y2 - y1;
            len2 = vx^2 + vy^2;

            if len2 == 0
                distToRay = sqrt((Xidx - x1).^2 + (Yidx - y1).^2);
            else
                t = ((Xidx - x1)*vx + (Yidx - y1)*vy) / len2;
                t = max(0, min(1, t));

                Xclosest = x1 + t*vx;
                Yclosest = y1 + t*vy;

                distToRay = sqrt((Xidx - Xclosest).^2 + ...
                                 (Yidx - Yclosest).^2);
            end

            mask4D(:,:,tx,rx) = distToRay <= halfWidth;
        end
    end
end

Nx_mask = size(mask4D, 1);
Ny_mask = size(mask4D, 2);

% Backward-compatible name for old plotting/debugging sections.
rayMask = mask4D;

% Example check:
% tx = 4; rx = 4;
% figure; plotXY(mask4D(:,:,tx,rx), xRange, yRange);
% title(sprintf('Mask tx %d to rx %d', tx, rx));

%% ------------------------------------------------------------------------

maxAngle = deg2rad(40);   % paper uses a maximum accepted angle around 40 deg

[epsr_final, epsr_num, epsr_den] = reconstructEpsrFinal( ...
    eps_r, angle, mask4D, maxAngle);

% Verification: if every TX-RX estimate is constant, overlaps must NOT
% amplify the image. The map should be constant wherever there is coverage.
% eps_r_test = ones(size(eps_r));
% [test_final, test_num, test_den] = reconstructEpsrFinal(eps_r_test, angle, mask4D, maxAngle);
% test_avg = averageEpsrFinal(test_num, test_den);
% disp([min(test_avg(:), [], 'omitnan'), max(test_avg(:), [], 'omitnan')]);

%% ------------------------------------------------------------------------

% Correct global weighted average over all transmitters and receivers.
% This is not the same as mean(epsr_final,3), because each pixel has a
% different number/weight of contributing masks.
epsr_avg = averageEpsrFinal(epsr_num, epsr_den);
% Pad epsr_avg back to full Nx x Ny size for visualization
epsr_avg_full = ones(Nx, Ny);   % background/PML shown as epsr = 1

epsr_avg_full(xRange, yRange) = epsr_avg;

otherRunTime = toc;
figure;
plotXY(epsr_avg_full, xRange, yRange);
title('Weighted average \epsilon_r reconstruction');
disp("Sim runtime: %d, other runtime: %d", simRunTime, otherRunTime);

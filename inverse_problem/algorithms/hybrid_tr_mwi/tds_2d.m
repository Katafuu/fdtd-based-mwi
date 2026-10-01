% tds_2d Run and plot the interactive hybrid TR/TDS workflow.
% A supplied cfg is reused; otherwise the default is built.

if ~exist('cfg', 'var')
    addpath(fileparts(mfilename('fullpath')), '-begin');
    cfg = build_cfg(struct());
end
assert(isstruct(cfg) && isscalar(cfg), ...
    'tds_2d:InvalidConfig', 'cfg must be a scalar struct.');

%% Circular antenna array plot
figure;
plotAntennaArray(cfg.antennas.pos, cfg.antennas.center, ...
    cfg.Nx, cfg.Ny, cfg.pml.thickness);
title('Circular antenna array');

% Assume itr_run(cfg) returns tr_result.focusMagFrame.
[tr_result, focusPoint] = tr_mwi(cfg); %#ok<NASGU>

%% Pair geometry: distance, global ray angle, and relative off-boresight angles
% distance(tx,rx): physical distance from antenna tx to antenna rx
% beta_tx(tx,rx): angle between tx boresight and outgoing ray tx->rx
% beta_rx(tx,rx): angle between rx boresight and incoming ray rx->tx
% pairWeight(tx,rx): cos(beta_tx)*cos(beta_rx), with nonpositive values removed
[distance, globalAngle, beta_tx, beta_rx, pairWeight] = buildCircularPairGeometry( ...
    cfg.antennas.pos, cfg.dx, cfg.dy, cfg.antennas.center);

maxAngle = deg2rad(40);
pairAccepted = abs(beta_tx) <= maxAngle & abs(beta_rx) <= maxAngle & isfinite(distance);
pairAccepted(eye(cfg.antennas.numAntennas) == 1) = false;

% Optional diagnostics
fprintf('Accepted TX-RX pairs: %d out of %d possible off-diagonal pairs.\n', ...
    nnz(pairAccepted), cfg.antennas.numAntennas*(cfg.antennas.numAntennas-1));

%% Run object/reference simulations and extract excess delay steps
delaysteps = nan(cfg.antennas.numAntennas, cfg.antennas.numAntennas);
arrival_obj_steps = nan(cfg.antennas.numAntennas, cfg.antennas.numAntennas);
arrival_inc_steps = nan(cfg.antennas.numAntennas, cfg.antennas.numAntennas);

for tx = 1:cfg.antennas.numAntennas
    cfg.source.samples(:) = 0;
    cfg.source.samples(tx, :) = ...
        cfg.source.func((0:cfg.Nt-1) .* cfg.dt);

    rxList = setdiff(1:cfg.antennas.numAntennas, tx);
    receiver_indices_this_tx = cfg.antennas.pos(rxList, :);

    % Object scan
    fdtdMexResult = fdtd_mex(cfg);

    % Incident/reference scan in background medium
    incident_cfg = fdtdmat.setBackgroundDefault(cfg, cfg.grid.background);
    fdtdMexResult_inc = fdtd_mex(incident_cfg);

    % Parse Ez into internal convention Ez(x,y,t)
    Ez = fdtdMexResult.Ez;

    Ez_inc = fdtdMexResult_inc.Ez;

    arrival_obj = getArrivalSteps(Ez, receiver_indices_this_tx, cfg);
    arrival_inc = getArrivalSteps(Ez_inc, receiver_indices_this_tx, incident_cfg);

    arrival_obj_steps(tx, rxList) = arrival_obj.';
    arrival_inc_steps(tx, rxList) = arrival_inc.';
    delaysteps(tx, rxList) = arrival_obj.' - arrival_inc.';

    fprintf('Finished transmitter %d / %d.\n', tx, cfg.antennas.numAntennas);
end

%% Compute pairwise average permittivity estimates
% delaysteps is already arrival_obj - arrival_inc in units of time steps.
% Therefore do NOT subtract 1 here.
Delta_t = delaysteps * cfg.dt;

clampNegativeDelays = true;
if clampNegativeDelays
    Delta_t(Delta_t < 0) = 0;
end

eps_r = (1 + (cfg.c0 .* Delta_t) ./ distance).^2;
eps_r(eye(cfg.antennas.numAntennas) == 1) = NaN;

%% Build paper-style footprint masks for the circular array
xRange = cfg.pml.thickness+1 : cfg.Nx-cfg.pml.thickness;
yRange = cfg.pml.thickness+1 : cfg.Ny-cfg.pml.thickness;

% Use footprint size based on circular antenna arc spacing as a reasonable
% first value. Tune Lfp_cells depending on desired coverage/resolution.
arcSpacingCells = 2*pi*cfg.antennas.radius / cfg.antennas.numAntennas;
Lfp_cells = max(1, round(arcSpacingCells/2));

focusBias = 1;

[mask4D, xRange, yRange] = buildFootprintMask( ...
    cfg.Nx, cfg.Ny, cfg.pml.thickness, cfg.antennas.pos, ...
    cfg.antennas.pos, Lfp_cells, focusPoint, focusBias, ...
    cfg.antennas.doiMask);

% Remove self-pair masks.
for ant = 1:cfg.antennas.numAntennas
    mask4D(:,:,ant,ant) = false;
end

fprintf('Built circular footprint masks with Lfp_cells = %d.\n', Lfp_cells);

% Plot the accepted TR-shifted footprints and their overlap count.
footprintFigure = figure('Name', 'Accepted TR-shifted footprint coverage', ...
    'NumberTitle', 'off', 'Visible', 'on');
plotCircularFootprintCoverage( ...
    mask4D, xRange, yRange, pairAccepted, cfg, focusPoint);
drawnow;

%% Reconstruct weighted epsr maps
[epsr_final, epsr_num, epsr_den] = reconstructEpsrFromWeights( ...
    eps_r, pairWeight, pairAccepted, mask4D);

epsr_avg = averageEpsrFinal(epsr_num, epsr_den);

%% Pad non-PML reconstruction back to full grid size for visualization
epsr_avg_full = cfg.grid.background.epsr .* ones(cfg.Nx, cfg.Ny);
epsr_avg_full(xRange, yRange) = epsr_avg;

reconstructionFigure = figure;
reconstructionAxes = axes('Parent', reconstructionFigure);
plotMwiReconstruction(reconstructionAxes, cfg, epsr_avg_full);

% Keep the footprint diagnostic visible when the script finishes.
figure(footprintFigure);
drawnow;

%% Optional verification: constant pair estimates should reconstruct constant values where covered
% eps_r_test = ones(size(eps_r));
% [~, test_num, test_den] = reconstructEpsrFromWeights(eps_r_test, pairWeight, pairAccepted, mask4D);
% test_avg = averageEpsrFinal(test_num, test_den);
% disp([min(test_avg(:), [], 'omitnan'), max(test_avg(:), [], 'omitnan')]);

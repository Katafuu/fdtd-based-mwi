% tds_2d_original Run one headless hybrid TR/TDS case and save its figures.
% Run build_cfg.m first, or use run_cases.m to construct cfg automatically.

assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'tds_2d_original:MissingConfig', ...
    'Run build_cfg.m before tds_2d_original.m.');
assert(exist('caseRunName', 'var') == 1 && ...
    (ischar(caseRunName) || (isstring(caseRunName) && isscalar(caseRunName))), ...
    'tds_2d_original:MissingRunName', ...
    'Set caseRunName (for example, ''run_0001'') before running this script.');

caseScriptDir = fileparts(mfilename('fullpath'));
caseOutputDir = fullfile(caseScriptDir, 'figs', char(caseRunName), 'original_tds');
if ~isfolder(caseOutputDir)
    mkdir(caseOutputDir);
end

%% Time reversal without nested plotting
trCfg = cfg;
trCfg.opts.enablePlot = false;
[tr_result, focusPoint] = tr_mwi(trCfg);

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

focusBias = 0;

[mask4D, xRange, yRange] = buildFootprintMask( ...
    cfg.Nx, cfg.Ny, cfg.pml.thickness, cfg.antennas.pos, ...
    cfg.antennas.pos, Lfp_cells, focusPoint, focusBias, ...
    cfg.antennas.doiMask);

% Remove self-pair masks.
for ant = 1:cfg.antennas.numAntennas
    mask4D(:,:,ant,ant) = false;
end

fprintf('Built circular footprint masks with Lfp_cells = %d.\n', Lfp_cells);


%% Reconstruct weighted epsr maps
[epsr_final, epsr_num, epsr_den] = reconstructEpsrFromWeights( ...
    eps_r, pairWeight, pairAccepted, mask4D);

epsr_avg = averageEpsrFinal(epsr_num, epsr_den);

%% Pad non-PML reconstruction back to full grid size for visualization
epsr_avg_full = cfg.grid.background.epsr .* ones(cfg.Nx, cfg.Ny);
epsr_avg_full(xRange, yRange) = epsr_avg;

plotData = struct();
plotData.trResult = tr_result;
plotData.focusPoint = focusPoint;
plotData.mask4D = mask4D;
plotData.xRange = xRange;
plotData.yRange = yRange;
plotData.pairAccepted = pairAccepted;
plotData.epsrAvgFull = epsr_avg_full;
plot_mwi(cfg, plotData, caseOutputDir);

%% Optional verification: constant pair estimates should reconstruct constant values where covered
% eps_r_test = ones(size(eps_r));
% [~, test_num, test_den] = reconstructEpsrFromWeights(eps_r_test, pairWeight, pairAccepted, mask4D);
% test_avg = averageEpsrFinal(test_num, test_den);
% disp([min(test_avg(:), [], 'omitnan'), max(test_avg(:), [], 'omitnan')]);

%% Iterative time reversal with the MEX TMz FDTD solver
% This script builds a centered-target test problem, runs forward scattering,
% and feeds each time-reversal receiver result into the next iteration.
clc;clear;
algorithmDir = fileparts(mfilename('fullpath'));
if isempty(algorithmDir)
    algorithmDir = pwd;
end
workspaceRoot = fileparts(fileparts(fileparts(algorithmDir)));
addpath(fullfile(workspaceRoot, 'buildLib'));
addpath(fullfile(algorithmDir, 'lib'));
addpath(fullfile(workspaceRoot, 'forward_solver', 'mex'));
if exist('fdtd_mex', 'file') ~= 3
    error('itr_simple:MissingFdtdMex', ...
        'Build forward_solver/mex/fdtd_mex before running this script.');
end

%% 0. Initialize cfg
cfg = struct();

%% 1. Basic grid parameters
cfg.c0 = 3e8;
cfg.mu0 = 4*pi*1e-7;
cfg.eps0 = 8.854e-12;
cfg.Nx = 400;
cfg.Ny = 400;
cfg.dx = 1e-3;
cfg.dy = 1e-3;
cfg.dt = 1 / (cfg.c0 * sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
cfg.Nt = 600;
cfg.sizeZ = 1;
cfg.snapshotStart = 0;
cfg.snapshotStride = 1;
cfg.returnEz = true;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = true;


%% 2. Grid
cfg.grid = struct();
cfg.grid.background = struct( ...
    'epsr', 1.0, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.0, ...
    'cond_m', 0.0);
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy], ...
    cfg.grid.background, [0 0]);

%% 3. PML
cfg.pml = fdtdpml.build_rectangularPML(cfg.grid, 20, 3, ...
    -(3 + 1) * log(1e-60) / ...
    (2 * sqrt(cfg.mu0/cfg.eps0) * 20 * cfg.dx), 1.0);
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = 1.0;
cfg.pml.ay = 1.0;
cfg.pml.az = 1.0;
cfg.pml.thickness = 60;
cfg.pml.m = 3;
cfg.pml.R = 1e-60;
cfg.pml.eta0 = sqrt(cfg.mu0/cfg.eps0);
cfg.pml.sigma_max = -(cfg.pml.m + 1) * log(cfg.pml.R) / ...
    (2 * cfg.pml.eta0 * cfg.pml.thickness * cfg.dx);
cfg.pml.sigma_e = max(cfg.pml.condx, cfg.pml.condy);

%% 4. Targets
cfg.targets = struct();
cfg.targets.name = 'circle';
cfg.targets.properties = struct( ...
    'center', [round(cfg.Nx/2), round(cfg.Ny/2)], ...
    'radius', 20);
cfg.targets.material = struct( ...
    'epsr', 5.0, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.1, ...
    'cond_m', 0.0);
cfg.targets.mask = fdtdgeom.shape_circle(cfg.grid, ...
    cfg.targets.properties.center, cfg.targets.properties.radius, 'index').mask;
cfg.grid = fdtdmat.applyRegion(cfg.grid, cfg.targets, cfg.targets.material);

%% 5. Antennas
cfg.antennas.numAntennas = 12;
cfg.antennas.txAntennas = 7; % Select the initial TR transmitting antennas.
cfg.antennas.pmlPadding = 5;
cfg.antennas.focusPadding = 20;

cfg.antennas.center = [round(cfg.Nx/2), round(cfg.Ny/2)];
cfg.antennas.radius = floor(min(cfg.grid.sizeXY)/2) - cfg.pml.thickness - cfg.antennas.pmlPadding;
cfg = buildCircularAntennaArrayIdx(cfg);


% cfg.antennas.orientation = 'top';
% cfg = buildPlanarAntennaIdx(cfg);

%% 6. Source
pulse_width = 20e-12;
t_delay = 4 * pulse_width;
amplitude = 10;
t = (0:cfg.Nt-1) .* cfg.dt;
source_normalization = max(abs(-((t - t_delay) ./ pulse_width.^2) .* ...
    exp(-0.5 .* ((t - t_delay) ./ pulse_width).^2)));
cfg.source.location = cfg.antennas.pos;
cfg.source.func = @(physicalTime) amplitude .* ...
    (-((physicalTime - t_delay) ./ pulse_width.^2) .* ...
    exp(-0.5 .* ((physicalTime - t_delay) ./ pulse_width).^2)) ./ source_normalization;

cfg.opts = struct();
cfg.opts.storeFieldHistory = true; % Retain sampled TR Ez history in tr_result.
cfg.filename = '';

%% 7. Iterative time reversal
incident_cfg = fdtdmat.setBackgroundDefault(cfg, cfg.grid.background);
if isfield(incident_cfg, 'targets')
    incident_cfg.targets = struct([]);
end
tx = cfg.source.func(t);
tx = tx(:).';

% Initialize the source matrix with the existing pulse at the active transmitter.
cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);
cfg.source.samples(cfg.antennas.txAntennas, :) = ...
    repmat(tx, numel(cfg.antennas.txAntennas), 1);
tx_signal = cfg.source.samples;
trace_opts = struct( ...
    'applyTemporalWindow', true, ...
    'temporalWindowTau', pulse_width, ...
    'normalizeTraces', true);

numIterations = 1;
focusFramesMag = zeros(cfg.Nx, cfg.Ny, numIterations);
focusFramesEntropy = zeros(cfg.Nx, cfg.Ny, numIterations);
% Iterate forward scattering, receiver injection, and feedback trace processing.
for iterNum = 1:numIterations
    fprintf('Iteration %d/%d: total-field MEX propagation...\n', ...
        iterNum, numIterations);
    cfg.source.samples = tx_signal;
    result_tot = fdtd_mex(cfg);

    fprintf('Iteration %d/%d: incident-field MEX propagation...\n', ...
        iterNum, numIterations);
    incident_cfg.source.samples = tx_signal;
    result_inc = fdtd_mex(incident_cfg);
    rx_scat = computeScatteredReceivers(result_tot, result_inc);
    [tx_tr, ~] = prepareTracesForTR(rx_scat, cfg.dt, trace_opts);

    fprintf('Iteration %d/%d: time-reversal MEX propagation...\n', ...
        iterNum, numIterations);
    tr_cfg = cfg;
    tr_cfg.source.samples = tx_tr;
    
    if (numIterations == 1) % collapses to original TR case
        tr_cfg.grid.cond_e = tr_cfg.grid.cond_e_bg;
        tr_cfg.grid.cond_m = tr_cfg.grid.cond_m_bg;
        tr_cfg.grid.epsr = tr_cfg.grid.epsr_bg;
    end

    trMexResult = fdtd_mex(tr_cfg, true);
    tr_result = analyzeTimeReversal(trMexResult, ...
        cfg.antennas.doiMask ~= 0, cfg.opts.storeFieldHistory, 1);
    focusFramesMag(:, :, iterNum) = tr_result.focus_mag_image;
    focusFramesEntropy(:, :, iterNum) = tr_result.focus_entropy_image;

    rx_after_tr = tr_result.rx_signals;
    [tr_data_rev, ~] = prepareTracesForTR(rx_after_tr, cfg.dt, trace_opts);
    cfg.source.samples = tr_data_rev;
    tx_signal = cfg.source.samples;
end

% Retain the final iterative focus frames for optional post-run plotting.
tr_result.focus_mag_image = focusFramesMag(:, :, end);
tr_result.focus_entropy_image = focusFramesEntropy(:, :, end);

function result = analyzeTimeReversal(mexResult, focusMask, storeHistory, historyStride)
%analyzeTimeReversal Compute focusing products from the complete MEX Ez history.

EzAll = mexResult.Ez;
numSteps = size(EzAll, 3);
imageMax = zeros(size(focusMask));
focusMagImage = zeros(size(focusMask));
focusMagEz = zeros(size(focusMask));
focusMagPeakValue = -inf;
focusMagStep = 1;
focusEntropyImage = zeros(size(focusMask));
focusEntropyEz = zeros(size(focusMask));
focusEntropyValue = inf;
focusEntropyStep = 1;
entropyAbsTol = 1e-15;
entropyRelTol = 1e-12;

for step = 1:numSteps
    EzFrame = EzAll(:, :, step);
    framePower = abs(EzFrame).^2;
    imageMax = max(imageMax, framePower);

    framePeakValue = max(framePower(focusMask));
    if framePeakValue > focusMagPeakValue
        focusMagPeakValue = framePeakValue;
        focusMagStep = step;
        focusMagEz = EzFrame;
        focusMagImage = framePower;
    end

    roiAmplitude = abs(EzFrame(focusMask));
    roiPeakAmplitude = max(roiAmplitude);
    if roiPeakAmplitude > entropyAbsTol
        entropyCutoff = max(entropyAbsTol, entropyRelTol .* roiPeakAmplitude);
        validAmplitude = roiAmplitude(roiAmplitude >= entropyCutoff);
        normalizedPower = (validAmplitude ./ roiPeakAmplitude).^2;
        entropyWeights = normalizedPower ./ sum(normalizedPower);
        frameEntropyValue = 1 ./ sum(entropyWeights.^2);
    else
        frameEntropyValue = inf;
    end
    if frameEntropyValue < focusEntropyValue
        focusEntropyValue = frameEntropyValue;
        focusEntropyStep = step;
        focusEntropyEz = EzFrame;
        focusEntropyImage = framePower;
    end
end

result = struct();
result.Ez_reversed = EzAll(:, :, end);
result.image_max = imageMax;
result.tr_mode = "receiverInjection";
result.loss_mode = "native";
result.num_steps = numSteps;
result.rx_signals = mexResult.rx_signals;
result.focus_mag_image = focusMagImage;
result.focus_mag_Ez = focusMagEz;
result.focus_mag_step = focusMagStep;
result.focus_mag_peak_value = focusMagPeakValue;
result.focus_entropy_image = focusEntropyImage;
result.focus_entropy_Ez = focusEntropyEz;
result.focus_entropy_step = focusEntropyStep;
result.focus_entropy_value = focusEntropyValue;
result.focus_mask = focusMask;

if storeHistory
    historySteps = unique([1:historyStride:numSteps, numSteps]);
    result.Ez_reversed_all = EzAll(:, :, historySteps);
    result.history_steps = historySteps;
    result.history_stride = historyStride;
    result.store_magnetic_history = false;
end
end

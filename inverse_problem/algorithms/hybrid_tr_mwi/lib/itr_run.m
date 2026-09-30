function tr_result = itr_run(cfg, SNR_dB)
%itr_run Run iterative time reversal using a caller-supplied FDTD cfg.
% Optional SNR_dB adds AWGN to the scattered receiver traces before TR.

addReceiverNoise = nargin >= 2;
if addReceiverNoise
    validateattributes(SNR_dB, {'numeric'}, ...
        {'real', 'scalar', 'finite'}, mfilename, 'SNR_dB');
    receiverNoisePowers = zeros(cfg.opts.numIterations, 1);
end

%% Repository / path setup
algorithmDir = fileparts(fileparts(mfilename('fullpath')));
if isempty(algorithmDir)
    algorithmDir = pwd;
end
algorithmsDir = fileparts(algorithmDir);
workspaceRoot = fileparts(fileparts(fileparts(algorithmDir)));

buildLibDir = fullfile(workspaceRoot, 'buildLib');
mexDir = fullfile(workspaceRoot, 'forward_solver', 'mex');
trLibDir = fullfile(algorithmsDir, 'time_reversal', 'lib');
mwiLibDir = fullfile(algorithmsDir, 'transmission_mwi', 'lib');
hybridLibDir = fullfile(algorithmDir, 'lib');

addpath(buildLibDir, '-end');
addpath(mexDir, '-end');
addpath(trLibDir, '-end');
addpath(mwiLibDir, '-end');
addpath(hybridLibDir, '-begin');
addpath(algorithmDir, '-begin');

if exist('fdtd_mex', 'file') ~= 3
    error('itr_run:MissingFdtdMex', ...
        'Build forward_solver/mex/fdtd_mex before running itr_run.');
end

%% Iterative time reversal
incident_cfg = fdtdmat.setBackgroundDefault(cfg, cfg.grid.background);
if isfield(incident_cfg, 'targets')
    incident_cfg.targets = struct([]);
end

% Validate the declared transmitters and caller-supplied source matrix.
if ~isfield(cfg.antennas, 'txAntennas') || isempty(cfg.antennas.txAntennas)
    error('itr_run:MissingTxAntennas', ...
        'cfg.antennas.txAntennas must select at least one transmitter.');
end
txAntennas = cfg.antennas.txAntennas;
if ~isnumeric(txAntennas) || ~isvector(txAntennas) || ...
        any(~isfinite(txAntennas)) || any(txAntennas ~= round(txAntennas)) || ...
        any(txAntennas < 1) || any(txAntennas > cfg.antennas.numAntennas)
    error('itr_run:InvalidTxAntennas', ...
        'cfg.antennas.txAntennas must contain valid antenna indices.');
end
if ~isfield(cfg, 'source') || ~isstruct(cfg.source) || ...
        ~isfield(cfg.source, 'samples') || ...
        ~isnumeric(cfg.source.samples) || ~isreal(cfg.source.samples) || ...
        issparse(cfg.source.samples) || ...
        ~isequal(size(cfg.source.samples), ...
        [cfg.antennas.numAntennas cfg.Nt]) || ...
        any(~isfinite(cfg.source.samples), 'all')
    error('itr_run:InvalidSourceSamples', ...
        ['cfg.source.samples must be a finite, full, real numeric ' ...
        'numAntennas-by-Nt matrix.']);
end
if ~any(cfg.source.samples ~= 0, 'all')
    error('itr_run:EmptySourceSamples', ...
        'cfg.source.samples must contain at least one active source row.');
end
tx_signal = double(cfg.source.samples);
trace_opts = struct( ...
    'applyTemporalWindow', cfg.opts.applyTemporalWindow, ...
    'temporalWindowTau', cfg.opts.temporalWindowTau, ...
    'normalizeTraces', cfg.opts.normalizeTraces);

focusFramesMag = zeros(cfg.Nx, cfg.Ny, cfg.opts.numIterations);
focusFramesEntropy = zeros(cfg.Nx, cfg.Ny, cfg.opts.numIterations);

% Iterate forward scattering, receiver injection, and feedback trace processing.
for iterNum = 1:cfg.opts.numIterations
    fprintf('Iteration %d/%d: total-field MEX propagation...\n', ...
        iterNum, cfg.opts.numIterations);
    cfg.source.samples = tx_signal;
    result_tot = fdtd_mex(cfg);

    fprintf('Iteration %d/%d: incident-field MEX propagation...\n', ...
        iterNum, cfg.opts.numIterations);
    incident_cfg.source.samples = tx_signal;
    result_inc = fdtd_mex(incident_cfg);
    rx_scat = computeScatteredReceivers(result_tot, result_inc);
    if addReceiverNoise
        signalPower = mean(rx_scat(:).^2);
        noisePower = signalPower * 10^(-SNR_dB/10);
        if ~isfinite(noisePower)
            error('itr_run:InvalidNoisePower', ...
                'The supplied SNR produces a nonfinite receiver noise power.');
        end
        receiverNoisePowers(iterNum) = noisePower;
        rx_scat = rx_scat + sqrt(noisePower) .* randn(size(rx_scat));
    end
    [tx_tr, ~] = prepareTracesForTR(rx_scat, cfg.dt, trace_opts);

    fprintf('Iteration %d/%d: time-reversal MEX propagation...\n', ...
        iterNum, cfg.opts.numIterations);
    tr_cfg = cfg;
    tr_cfg.source.samples = tx_tr;

    if (cfg.opts.numIterations == 1) % collapses to original TR case
        tr_cfg.grid.cond_e = tr_cfg.grid.cond_e_bg;
        tr_cfg.grid.cond_m = tr_cfg.grid.cond_m_bg;
        tr_cfg.grid.epsr = tr_cfg.grid.epsr_bg;
    end
    
    trMexResult = fdtd_mex(tr_cfg, true);
    tr_result = analyzeTimeReversal(trMexResult, ...
        cfg.antennas.doiMask ~= 0, ...
        cfg.opts.storeFieldHistory, cfg.opts.historyStride);
    focusFramesMag(:, :, iterNum) = tr_result.focus_mag_image;
    focusFramesEntropy(:, :, iterNum) = tr_result.focus_entropy_image;

    rx_after_tr = tr_result.rx_signals;
    [tr_data_rev, ~] = prepareTracesForTR(rx_after_tr, cfg.dt, trace_opts);
    cfg.source.samples = tr_data_rev;
    tx_signal = cfg.source.samples;
end

tr_result.focus_mag_image = focusFramesMag(:, :, end);
tr_result.focus_entropy_image = focusFramesEntropy(:, :, end);
tr_result.focusMagFrame = tr_result.focus_mag_Ez;
if addReceiverNoise
    tr_result.receiverNoisePowers = receiverNoisePowers;
end
end
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

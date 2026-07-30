% fbts Reconstruct lossless relative permittivity from Ez data by FBTS.
% Run build_cfg.m immediately before this script.

assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'fbts:MissingConfig', 'Run build_cfg.m before fbts.m.');
assert(exist('fdtd_mex', 'file') == 3, ...
    'fbts:MissingFdtdMex', ...
    'Build forward_solver/mex/fdtd_mex before running fbts.m.');

%% Paper source pulse and inversion setup
time = (0:cfg.Nt-1) .* cfg.dt;
tau = 0.125e-9;
cfg.source.func = @(t) ...
    (4 .* t.^3 ./ tau.^4 - t.^4 ./ tau.^5) .* exp(-t ./ tau);
sourcePulse = cfg.source.func(time);
T = time(end);
K = cos(pi .* time ./ (2 .* T));
K(end) = 0;

assert(isfield(cfg, 'deltaF') && isnumeric(cfg.deltaF) && ...
    isscalar(cfg.deltaF) && isfinite(cfg.deltaF) && cfg.deltaF > 0, ...
    'fbts:InvalidDeltaF', ...
    'cfg.deltaF must be a finite positive frequency resolution.');
sensitivityDownsampleFactor = 4;

numIterations = 15;
transmitters = cfg.antennas.txAntennas;
numTransmitters = numel(transmitters);
numReceivers = cfg.antennas.numAntennas;
doiMask = logical(cfg.antennas.doiMask);

epsrTrue = cfg.grid.epsr;
epsrBackground = cfg.grid.epsr_bg;
epsrEst = epsrBackground;
epsrLowerBound = 1.0;

[doiX, doiY] = find(doiMask);
xRange = min(doiX):max(doiX);
yRange = min(doiY):max(doiY);

%% Synthetic measurements from the true material
EzMeasured = zeros(numTransmitters, numReceivers, cfg.Nt);
fprintf("Generating synthetic measurements\n");
for transmitterIndex = 1:numTransmitters
    transmitter = transmitters(transmitterIndex);
    measurementCfg = cfg;
    measurementCfg.returnEz = false;
    measurementCfg.returnHx = false;
    measurementCfg.returnHy = false;
    measurementCfg.returnRxSignals = true;
    measurementCfg.source.samples(:) = 0;
    measurementCfg.source.samples(transmitter, :) = sourcePulse;
    measuredResult = fdtd_mex(measurementCfg);
    EzMeasured(transmitterIndex, :, :) = reshape( ...
        measuredResult.rx_signals, 1, numReceivers, cfg.Nt);

    clear measuredResult
    fprintf('Measured transmitter %d / %d.\n', ...
        transmitterIndex, numTransmitters);
end

%% Forward-backward inversion
costHistory = zeros(numIterations, 1);
iterationRuntime = zeros(numIterations, 1);
relativeErrorPercent = zeros(numIterations, 1);
alphaHistory = zeros(numIterations, 1);
prBetaHistory = zeros(numIterations, 1);
prRestartHistory = false(numIterations, 1);
projectionAppliedHistory = false(numIterations, 1);
projectedCellCountHistory = zeros(numIterations, 1);
directionalDerivativeHistory = zeros(numIterations, 1);
stepAHistory = zeros(numIterations, 1);
stepQHistory = zeros(numIterations, 1);
finiteDifferenceHHistory = zeros(numIterations, 1);
sensitivityDtHistory = zeros(numIterations, 1);
sensitivityNtHistory = zeros(numIterations, 1);
EzModel = zeros(numTransmitters, numReceivers, cfg.Nt);
EzSensitivity = zeros(numTransmitters, numReceivers, cfg.Nt);
previousProjectedGradEpsr = zeros(size(epsrEst));
previousDirectionEpsr = zeros(size(epsrEst));
directionEpsr = zeros(size(epsrEst));

fprintf("Starting FBTS loop\n");
for iteration = 1:numIterations
    iterationTimer = tic;
    gradEpsr = zeros(size(epsrEst));
    iterationCost = 0;

    for transmitterIndex = 1:numTransmitters
        transmitter = transmitters(transmitterIndex);

        %% Forward run in the current estimate
        forwardCfg = cfg;
        forwardCfg.grid.epsr = epsrEst;
        forwardCfg.returnEz = true;
        forwardCfg.returnHx = false;
        forwardCfg.returnHy = false;
        forwardCfg.returnRxSignals = true;
        forwardCfg.source.samples(:) = 0;
        forwardCfg.source.samples(transmitter, :) = sourcePulse;
        forwardResult = fdtd_mex(forwardCfg);

        modelEz = forwardResult.rx_signals;
        EzModel(transmitterIndex, :, :) = reshape( ...
            modelEz, 1, numReceivers, cfg.Nt);
        EzForward = forwardResult.Ez(xRange, yRange, :);
        clear forwardResult

        %% Weighted model-minus-measured Ez residual
        residualEz = modelEz - ...
            reshape(EzMeasured(transmitterIndex, :, :), ...
            numReceivers, cfg.Nt);
        weightedResidualEz = K .* residualEz;
        iterationCost = iterationCost + ...
            sum(K .* residualEz.^2, 'all') .* cfg.dt;

        %% Native C time-reversal run with all receiver residuals
        adjointCfg = forwardCfg;
        adjointCfg.returnRxSignals = false;
        adjointCfg.source.samples = -flip(weightedResidualEz, 2);
        adjointResult = fdtd_mex(adjointCfg, true);
        wHatM1TR = adjointResult.Ez(xRange, yRange, :);
        clear adjointResult

        %% Restore forward physical-time order
        wHatM1 = flip(wHatM1TR, 3);

        %% Permittivity gradient, equation (32)
        dEzDct = diff(EzForward, 1, 3) ./ (cfg.c0 .* cfg.dt);
        w1Half = 0.5 .* ...
            (wHatM1(:, :, 2:end) + wHatM1(:, :, 1:end-1));
        gradEpsrM = 2 .* sum(w1Half .* dEzDct, 3) .* cfg.dt;
        gradEpsr(xRange, yRange) = ...
            gradEpsr(xRange, yRange) + gradEpsrM;

        clear EzForward wHatM1TR wHatM1
        fprintf('Iteration %d/%d, transmitter %d/%d.\n', ...
            iteration, numIterations, transmitterIndex, numTransmitters);
    end

    %% Bound-projected Polak-Ribiere-Polyak direction, equation (23)
    gradEpsr(~doiMask) = 0;
    projectedGradEpsr = gradEpsr;
    activeLowerBound = doiMask & epsrEst <= epsrLowerBound;
    projectedGradEpsr(activeLowerBound & projectedGradEpsr > 0) = 0;

    restartFromProjection = iteration > 1 && ...
        projectionAppliedHistory(iteration - 1);
    if iteration == 1 || restartFromProjection
        prBeta = 0;
        directionEpsr = -projectedGradEpsr;
        prRestart = restartFromProjection;
    else
        previousGradientNormSquared = ...
            sum(previousProjectedGradEpsr(doiMask).^2, 'all');
        assert(isfinite(previousGradientNormSquared) && ...
            previousGradientNormSquared > 0, ...
            'fbts:InvalidPreviousGradient', ...
            'The previous DOI gradient must have a finite positive norm.');

        gradientChange = projectedGradEpsr - previousProjectedGradEpsr;
        prBeta = sum(gradientChange(doiMask) .* ...
            projectedGradEpsr(doiMask), 'all') ./ ...
            previousGradientNormSquared;
        directionEpsr = -projectedGradEpsr + ...
            prBeta .* previousDirectionEpsr;
        prRestart = false;
    end
    directionEpsr(~doiMask) = 0;

    directionalDerivative = sum(projectedGradEpsr(doiMask) .* ...
        directionEpsr(doiMask), 'all');
    directionNeedsRestart = ~isfinite(prBeta) || ...
        any(~isfinite(directionEpsr(doiMask)), 'all') || ...
        ~isfinite(directionalDerivative) || directionalDerivative >= 0;
    if directionNeedsRestart
        prBeta = 0;
        directionEpsr = -projectedGradEpsr;
        directionEpsr(~doiMask) = 0;
        directionalDerivative = sum(projectedGradEpsr(doiMask) .* ...
            directionEpsr(doiMask), 'all');
        prRestart = true;
    end

    assert(isfinite(directionalDerivative) && directionalDerivative < 0 && ...
        all(isfinite(directionEpsr(doiMask)), 'all'), ...
        'fbts:InvalidSearchDirection', ...
        'The projected search direction must be finite and descending.');

    %% Directional finite difference for the Frechet field, equation (24)
    directionMaximum = max(abs(directionEpsr(doiMask)), [], 'all');
    assert(isfinite(directionMaximum) && directionMaximum > 0, ...
        'fbts:ZeroSearchDirection', ...
        'The DOI search direction must have a finite nonzero magnitude.');
    finiteDifferenceH = 1e-3 ./ directionMaximum;
    assert(all(isfinite(finiteDifferenceH) & finiteDifferenceH > 0, 'all'), ...
        'fbts:InvalidFiniteDifferenceStep', ...
        'The feasible finite-difference step must be finite and positive.');
    epsrPerturbed = epsrEst + finiteDifferenceH .* directionEpsr;
    epsrPerturbed(doiMask) = max( ...
        epsrPerturbed(doiMask), epsrLowerBound);
    epsrPerturbed(~doiMask) = epsrBackground(~doiMask);

    stepA = 0;
    stepQ = 0;

    coarseBaseCfg = cfg;
    coarseBaseCfg.grid.epsr = epsrEst;
    coarseBaseCfg = prepareCoarseSensitivityCfg( ...
        coarseBaseCfg, sensitivityDownsampleFactor);
    coarsePerturbedCfg = cfg;
    coarsePerturbedCfg.grid.epsr = epsrPerturbed;
    coarsePerturbedCfg = prepareCoarseSensitivityCfg( ...
        coarsePerturbedCfg, sensitivityDownsampleFactor);
    assert(coarseBaseCfg.dt == coarsePerturbedCfg.dt && ...
        coarseBaseCfg.Nt == coarsePerturbedCfg.Nt, ...
        'fbts:SensitivityTimeGridMismatch', ...
        'Coarse baseline and perturbation time grids must match.');
    coarseTime = (0:coarseBaseCfg.Nt-1) .* coarseBaseCfg.dt;
    coarseSourcePulse = cfg.source.func(coarseTime);
    coarseK = interp1(time, K, coarseTime, 'linear', 0);
    sensitivityDtHistory(iteration) = coarseBaseCfg.dt;
    sensitivityNtHistory(iteration) = coarseBaseCfg.Nt;

    fprintf('Computing iteration %d directional sensitivity\n', iteration);
    for transmitterIndex = 1:numTransmitters
        transmitter = transmitters(transmitterIndex);
        modelEz = reshape(EzModel(transmitterIndex, :, :), ...
            numReceivers, cfg.Nt);
        residualEz = modelEz - ...
            reshape(EzMeasured(transmitterIndex, :, :), ...
            numReceivers, cfg.Nt);
        coarseResidualEz = interp1( ...
            time, residualEz.', coarseTime, 'linear', 0).';

        coarseBaseCfg.source.samples(:) = 0;
        coarseBaseCfg.source.samples(transmitter, :) = coarseSourcePulse;
        coarseBaselineResult = fdtd_mex(coarseBaseCfg);
        coarsePerturbedCfg.source.samples(:) = 0;
        coarsePerturbedCfg.source.samples(transmitter, :) = coarseSourcePulse;
        coarsePerturbedResult = fdtd_mex(coarsePerturbedCfg);

        coarseSensitivityEz = ...
            (coarsePerturbedResult.rx_signals - ...
            coarseBaselineResult.rx_signals) ./ finiteDifferenceH;
        sensitivityEz = interp1( ...
            coarseTime, coarseSensitivityEz.', time, 'linear', 0).';
        EzSensitivity(transmitterIndex, :, :) = reshape( ...
            sensitivityEz, 1, numReceivers, cfg.Nt);
        stepA = stepA + ...
            sum(coarseK .* coarseSensitivityEz.^2, 'all') .* ...
            coarseBaseCfg.dt;
        stepQ = stepQ - ...
            sum(coarseK .* coarseResidualEz .* coarseSensitivityEz, 'all') .* ...
            coarseBaseCfg.dt;

        clear coarseBaselineResult coarsePerturbedResult
        fprintf('Sensitivity iteration %d/%d, transmitter %d/%d.\n', ...
            iteration, numIterations, transmitterIndex, numTransmitters);
    end

    %% Data-optimal step size, equations (31)-(33)
    assert(isfinite(stepA) && stepA > 0 && isfinite(stepQ), ...
        'fbts:InvalidOptimalStepSystem', ...
        'The optimal-step terms require finite q and finite positive a.');
    unconstrainedAlpha = stepQ ./ stepA;
    assert(isfinite(unconstrainedAlpha), 'fbts:InvalidOptimalStep', ...
        'The optimal permittivity step must be finite.');
    alpha = max(unconstrainedAlpha, 0);

    updatedEpsrDoi = epsrEst(doiMask) + ...
        alpha .* directionEpsr(doiMask);
    projectedCellCount = nnz(updatedEpsrDoi < epsrLowerBound);
    projectionApplied = projectedCellCount > 0;
    epsrEst(doiMask) = max(updatedEpsrDoi, epsrLowerBound);
    epsrEst(~doiMask) = epsrBackground(~doiMask);

    previousProjectedGradEpsr = projectedGradEpsr;
    previousDirectionEpsr = directionEpsr;
    alphaHistory(iteration) = alpha;
    prBetaHistory(iteration) = prBeta;
    prRestartHistory(iteration) = prRestart;
    projectionAppliedHistory(iteration) = projectionApplied;
    projectedCellCountHistory(iteration) = projectedCellCount;
    directionalDerivativeHistory(iteration) = directionalDerivative;
    stepAHistory(iteration) = stepA;
    stepQHistory(iteration) = stepQ;
    finiteDifferenceHHistory(iteration) = finiteDifferenceH;


    relativeErrorPercent(iteration) = 100 .* ...
        norm(epsrEst(doiMask) - epsrTrue(doiMask)) ./ ...
        norm(epsrTrue(doiMask));

    costHistory(iteration) = iterationCost;
    fprintf('Iteration %d cost: %.6e\n', iteration, iterationCost);
    fprintf(['Iteration %d optimal step: alpha = %.6e, ' ...
        'beta = %.6e, a = %.6e, q = %.6e, h = %.6e.\n'], ...
        iteration, alpha, prBeta, stepA, stepQ, finiteDifferenceH);
    fprintf('Iteration %d/%d relative DOI error: %.6e%%.\n', ...
        iteration, numIterations, relativeErrorPercent(iteration));
    iterationRuntime(iteration) = toc(iterationTimer);
    fprintf('Iteration %d/%d runtime: %.3f seconds.\n', ...
        iteration, numIterations, iterationRuntime(iteration));
end

totalRuntime = sum(iterationRuntime);
averageRuntime = mean(iterationRuntime);
fprintf('Total FBTS runtime: %.3f seconds.\n', totalRuntime);
fprintf('Average FBTS iteration runtime: %.3f seconds.\n', averageRuntime);

%% Results
results = struct();
results.epsr_true = epsrTrue;
results.epsr_est = epsrEst;
results.gradient_epsr = projectedGradEpsr;
results.raw_gradient_epsr = gradEpsr;
results.direction_epsr = directionEpsr;
results.cost = costHistory;
results.alpha = alphaHistory;
results.pr_beta = prBetaHistory;
results.pr_restart = prRestartHistory;
results.projection_applied = projectionAppliedHistory;
results.projected_cell_count = projectedCellCountHistory;
results.directional_derivative = directionalDerivativeHistory;
results.epsr_lower_bound = epsrLowerBound;
results.step_a = stepAHistory;
results.step_q = stepQHistory;
results.finite_difference_h = finiteDifferenceHHistory;
results.sensitivity_downsample_factor = sensitivityDownsampleFactor;
results.sensitivity_dt = sensitivityDtHistory;
results.sensitivity_Nt = sensitivityNtHistory;
results.iteration_runtime = iterationRuntime;
results.total_runtime = totalRuntime;
results.average_runtime = averageRuntime;
results.iteration = (1:numIterations).';
results.relative_error_percent = relativeErrorPercent;
results.relative_error_definition = [ ...
    '100 * norm(epsr_est(doiMask) - epsr_true(doiMask)) / ' ...
    'norm(epsr_true(doiMask))'];
results.doi_mask = doiMask;
results.Ez_measured = EzMeasured;
results.Ez_model = EzModel;
results.Ez_sensitivity = EzSensitivity;

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end
runOutputDirectory = createNextRunDirectory(fullfile(scriptDir, 'figs'));
figureFiles = [ ...
    string(fullfile(runOutputDirectory, ...
        '01_true_relative_permittivity.png'))
    string(fullfile(runOutputDirectory, ...
        '02_estimated_relative_permittivity.png'))
    string(fullfile(runOutputDirectory, ...
        '03_relative_error_convergence.png'))
    string(fullfile(runOutputDirectory, ...
        '04_cost_convergence.png'))
];
resultFile = string(fullfile(runOutputDirectory, 'fbts_run_data.mat'));
results.output_directory = string(runOutputDirectory);
results.output_files = struct( ...
    'figures', figureFiles, ...
    'mat_file', resultFile);

plot_fbts;
save(char(resultFile), 'cfg', 'results', '-v7.3');
fprintf('Saved FBTS figures and run data under %s.\n', ...
    runOutputDirectory);

function runDirectory = createNextRunDirectory(outputRoot)
outputRoot = char(outputRoot);
if ~isfolder(outputRoot)
    [created, message] = mkdir(outputRoot);
    if ~created
        error('fbts:CreateOutputRootFailed', '%s', message);
    end
end

existingRuns = dir(fullfile(outputRoot, 'run_*'));
existingRuns = existingRuns([existingRuns.isdir]);
runNumbers = nan(numel(existingRuns), 1);
numMatchedRuns = 0;
for directoryIndex = 1:numel(existingRuns)
    token = regexp(existingRuns(directoryIndex).name, ...
        '^run_(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        numMatchedRuns = numMatchedRuns + 1;
        runNumbers(numMatchedRuns) = str2double(token{1});
    end
end
runNumbers = runNumbers(1:numMatchedRuns);
if isempty(runNumbers)
    nextRunNumber = 1;
else
    nextRunNumber = max(runNumbers) + 1;
end

runDirectory = fullfile(outputRoot, sprintf('run_%04d', nextRunNumber));
[created, message] = mkdir(runDirectory);
if ~created
    error('fbts:CreateRunDirectoryFailed', '%s', message);
end
end

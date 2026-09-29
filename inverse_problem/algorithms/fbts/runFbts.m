function [results, cfg] = runFbts(cfg, fbtsOptions, measurements, execution)
%runFbts Reconstruct lossless relative permittivity without plotting or saving.
%   [results, cfg] = runFbts(cfg) uses the defaults in resolveFbtsOptions.
%   Optional fbtsOptions controls iteration count and sensitivity coarsening.
%   The returned cfg includes the source pulse used for the reconstruction.
%   Result costs/signals precede each update; DOI errors follow each update.

setupTimer = tic;
assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'fbts:MissingConfig', 'Run build_cfg.m before fbts_demo.m.');
assert(exist('fdtd_mex', 'file') == 3, ...
    'fbts:MissingFdtdMex', ...
    'Build forward_solver/mex/fdtd_mex before running fbts_demo.m.');

if nargin < 2
    fbtsOptions = [];
end
activeFbtsOptions = resolveFbtsOptions(fbtsOptions);
cfg.solverThreads = activeFbtsOptions.solverThreads;
cfg.solverBackend = activeFbtsOptions.solverBackend;
cfg.gpuDeviceIndex = activeFbtsOptions.gpuDeviceIndex;
if nargin < 4, execution = []; end
execution = resolveFbtsExecution(execution);

%% Paper source pulse and inversion setup
[cfg, time, sourcePulse, K] = prepareFbtsSource(cfg);

assert(isfield(cfg, 'deltaF') && isnumeric(cfg.deltaF) && ...
    isscalar(cfg.deltaF) && isfinite(cfg.deltaF) && cfg.deltaF > 0, ...
    'fbts:InvalidDeltaF', ...
    'cfg.deltaF must be a finite positive frequency resolution.');
sensitivityDownsampleFactor = activeFbtsOptions.sensitivityDownsampleFactor;

numIterations = activeFbtsOptions.numIterations;
transmitters = cfg.antennas.txAntennas;
numTransmitters = numel(transmitters);
numReceivers = cfg.antennas.numAntennas;
doiMask = logical(cfg.antennas.doiMask);

% Rows identify transmitters; columns identify FBTS iterations. Column-major
% order is therefore chronological when collecting per-phase statistics.
fdtdTiming = struct( ...
    'measurement', zeros(numTransmitters, 1), ...
    'forward', nan(numTransmitters, numIterations), ...
    'adjoint', nan(numTransmitters, numIterations), ...
    'sensitivity_baseline', nan(numTransmitters, numIterations), ...
    'sensitivity_perturbation', nan(numTransmitters, numIterations));

phaseTiming = struct();
for name = {'residual','gradient','sensitivity_processing'}
    phaseTiming.(name{1}) = nan(numTransmitters,numIterations);
end
for name = {'direction','coarse_setup','update','history_write','checkpoint_write','history_hash','iteration_wall'}
    phaseTiming.(name{1}) = nan(numIterations,1);
end
s = struct();
epsrTrue = cfg.grid.epsr;
epsrBackground = cfg.grid.epsr_bg;
s.epsrEst = epsrBackground;
epsrLowerBound = 1.0;

[doiX, doiY] = find(doiMask);
xRange = min(doiX):max(doiX);
yRange = min(doiY):max(doiY);

%% Measurements may come from the full-array scan of this same scene.
if nargin < 3 || isempty(measurements)
    measurements = acquireFbtsMeasurements(cfg);
    fdtdTiming.measurement = measurements.timing.solver_seconds;
else
    fdtdTiming.measurement = zeros(0,1);
end
EzMeasured = selectFbtsMeasurements(measurements, cfg);

%% Forward-backward inversion
h = struct();
h.cost = zeros(numIterations, 1);
h.iteration_runtime = zeros(numIterations, 1);
h.relative_error_percent = zeros(numIterations, 1);
h.alpha = zeros(numIterations, 1);
h.pr_beta = zeros(numIterations, 1);
h.pr_restart = false(numIterations, 1);
h.projection_applied = false(numIterations, 1);
h.projected_cell_count = zeros(numIterations, 1);
h.directional_derivative = zeros(numIterations, 1);
h.step_a = zeros(numIterations, 1);
h.step_q = zeros(numIterations, 1);
h.finite_difference_h = zeros(numIterations, 1);
h.sensitivity_dt = zeros(numIterations, 1);
h.sensitivity_Nt = zeros(numIterations, 1);
EzModel = zeros(numTransmitters, numReceivers, cfg.Nt);
EzSensitivity = zeros(numTransmitters, numReceivers, cfg.Nt);
s.previousProjectedGradEpsr = zeros(size(s.epsrEst));
s.previousDirectionEpsr = zeros(size(s.epsrEst));
s.directionEpsr = zeros(size(s.epsrEst));

h.cost_by_tx = zeros(numTransmitters,numIterations);
h.unconstrained_alpha = zeros(numIterations,1);
h.raw_gradient_norm = zeros(numIterations,1);
h.projected_gradient_norm = zeros(numIterations,1);
h.update_norm = zeros(numIterations,1);
h.restart_reason = repmat("",numIterations,1);
s.gradEpsr = zeros(size(s.epsrEst));
s.projectedGradEpsr = zeros(size(s.epsrEst));
s.lastModel = EzModel; s.lastSensitivity = EzSensitivity;
coarseTemplate = cfg;
coarseTemplate.grid.epsr = epsrBackground;
coarseTemplate = prepareCoarseSensitivityCfg(coarseTemplate,sensitivityDownsampleFactor);
provenance = fbtsCodeProvenance();
assert(isempty(execution.expectedCodeHash) || strcmp(execution.expectedCodeHash,provenance.code_hash), ...
    'fbts:CodeChanged','Source/MEX code changed after this batch was prepared.');
runIdentity = fbtsHash({fbtsPortableConfig(cfg), numIterations, ...
    sensitivityDownsampleFactor, measurements.metadata.fingerprint, ...
    provenance.code_hash, execution.computeImageError});
firstIteration = 1;
if execution.resume && ~isempty(execution.checkpointFile) && isfile(execution.checkpointFile)
    loaded = load(execution.checkpointFile,'checkpoint');
    cp = loaded.checkpoint;
    assert(strcmp(cp.checksum,fbtsHash(rmfield(cp,'checksum'))) && ...
        strcmp(cp.identity,runIdentity), 'fbts:CheckpointMismatch', ...
        'Checkpoint inputs, code, or checksum differ.');
    s = cp.state; h = cp.history; fdtdTiming = cp.fdtdTiming;
    if isfield(cp,'phaseTiming'), phaseTiming = cp.phaseTiming; end
    if ~isfield(h,'relative_error_percent'), h.relative_error_percent = NaN(numIterations,1); end
    % Reusing a scan on resume does not repeat its acquisition work.
    fdtdTiming.measurement = zeros(0,1);
    firstIteration = cp.completed_iterations+1;
    assert(isfile(execution.historyFile), 'fbts:MissingHistory', 'Checkpoint history is missing.');
    traceHistory = matfile(execution.historyFile,'Writable',true);
    assert(strcmp(traceHistory.identity,runIdentity) && ...
        traceHistory.complete_iterations >= cp.completed_iterations, ...
        'fbts:CheckpointMismatch','Checkpoint and trace history are inconsistent.');
    assert(isequal(traceHistory.epsr_history_doi(:,firstIteration),s.epsrEst(doiMask)), ...
        'fbts:CheckpointMismatch','Checkpoint image differs from history.');
    assert(strcmp(cp.trace_checksum,historyChecksum(traceHistory,firstIteration-1,numIterations)), ...
        'fbts:CheckpointMismatch','Committed iteration history is modified.');
else
    traceHistory = struct('identity',runIdentity,'complete_iterations',0, ...
        'epsr_history_doi',zeros(nnz(doiMask),numIterations+1), ...
        'rx_model_history',zeros(numTransmitters,numReceivers,cfg.Nt,numIterations), ...
        'rx_sensitivity_coarse_history',zeros(numTransmitters,numReceivers,coarseTemplate.Nt,numIterations));
    traceHistory.epsr_history_doi(:,1) = s.epsrEst(doiMask);
    if ~isempty(execution.historyFile)
        saveFbtsAtomic(execution.historyFile,traceHistory);
        traceHistory = matfile(execution.historyFile,'Writable',true);
    end
end

phaseTiming.setup_seconds = toc(setupTimer);
completedIterations = firstIteration-1;
iteration = firstIteration; transmitterIndex = NaN; phase = 'iteration_setup';
fprintf("Starting FBTS loop\n");
try
for iteration = firstIteration:numIterations
    iterationTimer = tic;
    s.gradEpsr = zeros(size(s.epsrEst));
    iterationCost = 0;
    previousEstimateDoi = s.epsrEst(doiMask);
    coarseSensitivityAll = zeros(numTransmitters,numReceivers,coarseTemplate.Nt);

    for transmitterIndex = 1:numTransmitters
        transmitter = transmitters(transmitterIndex);

        %% Forward run in the current estimate
        forwardCfg = cfg;
        forwardCfg.grid.epsr = s.epsrEst;
        forwardCfg.returnEz = true;
        forwardCfg.returnHx = false;
        forwardCfg.returnHy = false;
        forwardCfg.returnRxSignals = true;
        forwardCfg.ezRegion = [xRange(1),xRange(end),yRange(1),yRange(end)];
        forwardCfg.source.samples(:) = 0;
        forwardCfg.source.samples(transmitter, :) = sourcePulse;
        phase = 'forward';
        fdtdTimer = tic;
        forwardResult = fbtsSolve(forwardCfg);
        fdtdTiming.forward(transmitterIndex, iteration) = toc(fdtdTimer);

        modelEz = forwardResult.rx_signals;
        EzModel(transmitterIndex, :, :) = reshape( ...
            modelEz, 1, numReceivers, cfg.Nt);
        EzForward = forwardResult.Ez;
        clear forwardResult

        phaseTimer = tic;
        %% Weighted model-minus-measured Ez residual
        residualEz = modelEz - ...
            reshape(EzMeasured(transmitterIndex, :, :), ...
            numReceivers, cfg.Nt);
        weightedResidualEz = K .* residualEz;
        txCost = sum(K .* residualEz.^2, 'all') .* cfg.dt;
        h.cost_by_tx(transmitterIndex,iteration) = txCost;
        iterationCost = iterationCost + txCost;

        phaseTiming.residual(transmitterIndex,iteration) = toc(phaseTimer);
        %% Native C time-reversal run with all receiver residuals
        adjointCfg = forwardCfg;
        adjointCfg.returnRxSignals = false;
        adjointCfg.source.samples = -flip(weightedResidualEz, 2);
        phase = 'adjoint';
        fdtdTimer = tic;
        adjointResult = fbtsSolve(adjointCfg, true);
        fdtdTiming.adjoint(transmitterIndex, iteration) = toc(fdtdTimer);
        wHatM1TR = adjointResult.Ez;
        clear adjointResult

        phaseTimer = tic;
        %% Restore forward physical-time order
        wHatM1 = flip(wHatM1TR, 3);

        %% Permittivity gradient, equation (32)
        dEzDct = diff(EzForward, 1, 3) ./ (cfg.c0 .* cfg.dt);
        w1Half = 0.5 .* ...
            (wHatM1(:, :, 2:end) + wHatM1(:, :, 1:end-1));
        gradEpsrM = 2 .* sum(w1Half .* dEzDct, 3) .* cfg.dt;
        s.gradEpsr(xRange, yRange) = ...
            s.gradEpsr(xRange, yRange) + gradEpsrM;

        phaseTiming.gradient(transmitterIndex,iteration) = toc(phaseTimer);
        clear EzForward wHatM1TR wHatM1 dEzDct w1Half gradEpsrM
        fprintf('Iteration %d/%d, transmitter %d/%d.\n', ...
            iteration, numIterations, transmitterIndex, numTransmitters);
    end

    phase = 'gradient_and_direction'; phaseTimer = tic;
    %% Bound-projected Polak-Ribiere-Polyak direction, equation (23)
    s.gradEpsr(~doiMask) = 0;
    s.projectedGradEpsr = s.gradEpsr;
    activeLowerBound = doiMask & s.epsrEst <= epsrLowerBound;
    s.projectedGradEpsr(activeLowerBound & s.projectedGradEpsr > 0) = 0;

    h.raw_gradient_norm(iteration) = norm(s.gradEpsr(doiMask));
    h.projected_gradient_norm(iteration) = norm(s.projectedGradEpsr(doiMask));
    restartFromProjection = iteration > 1 && ...
        h.projection_applied(iteration - 1);
    if iteration == 1 || restartFromProjection
        prBeta = 0;
        s.directionEpsr = -s.projectedGradEpsr;
        prRestart = restartFromProjection;
        if restartFromProjection, h.restart_reason(iteration) = "projection"; end
    else
        previousGradientNormSquared = ...
            sum(s.previousProjectedGradEpsr(doiMask).^2, 'all');
        assert(isfinite(previousGradientNormSquared) && ...
            previousGradientNormSquared > 0, ...
            'fbts:InvalidPreviousGradient', ...
            'The previous DOI gradient must have a finite positive norm.');

        gradientChange = s.projectedGradEpsr - s.previousProjectedGradEpsr;
        prBeta = sum(gradientChange(doiMask) .* ...
            s.projectedGradEpsr(doiMask), 'all') ./ ...
            previousGradientNormSquared;
        s.directionEpsr = -s.projectedGradEpsr + ...
            prBeta .* s.previousDirectionEpsr;
        prRestart = false;
    end
    s.directionEpsr(~doiMask) = 0;

    directionalDerivative = sum(s.projectedGradEpsr(doiMask) .* ...
        s.directionEpsr(doiMask), 'all');
    directionNeedsRestart = ~isfinite(prBeta) || ...
        any(~isfinite(s.directionEpsr(doiMask)), 'all') || ...
        ~isfinite(directionalDerivative) || directionalDerivative >= 0;
    if directionNeedsRestart
        prBeta = 0;
        s.directionEpsr = -s.projectedGradEpsr;
        s.directionEpsr(~doiMask) = 0;
        directionalDerivative = sum(s.projectedGradEpsr(doiMask) .* ...
            s.directionEpsr(doiMask), 'all');
        prRestart = true;
        h.restart_reason(iteration) = "non_descent_or_nonfinite";
    end

    assert(isfinite(directionalDerivative) && directionalDerivative < 0 && ...
        all(isfinite(s.directionEpsr(doiMask)), 'all'), ...
        'fbts:InvalidSearchDirection', ...
        'The projected search direction must be finite and descending.');

    %% Directional finite difference for the Frechet field, equation (24)
    directionMaximum = max(abs(s.directionEpsr(doiMask)), [], 'all');
    assert(isfinite(directionMaximum) && directionMaximum > 0, ...
        'fbts:ZeroSearchDirection', ...
        'The DOI search direction must have a finite nonzero magnitude.');
    finiteDifferenceH = 1e-3 ./ directionMaximum;
    assert(all(isfinite(finiteDifferenceH) & finiteDifferenceH > 0, 'all'), ...
        'fbts:InvalidFiniteDifferenceStep', ...
        'The feasible finite-difference step must be finite and positive.');
    epsrPerturbed = s.epsrEst + finiteDifferenceH .* s.directionEpsr;
    epsrPerturbed(doiMask) = max( ...
        epsrPerturbed(doiMask), epsrLowerBound);
    epsrPerturbed(~doiMask) = epsrBackground(~doiMask);

    stepA = 0;
    stepQ = 0;

    phaseTiming.direction(iteration) = toc(phaseTimer);
    phaseTimer = tic;
    coarseBaseCfg = cfg;
    coarseBaseCfg.grid.epsr = s.epsrEst;
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
    h.sensitivity_dt(iteration) = coarseBaseCfg.dt;
    h.sensitivity_Nt(iteration) = coarseBaseCfg.Nt;

    phaseTiming.coarse_setup(iteration) = toc(phaseTimer);
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
        phase = 'sensitivity_baseline';
        fdtdTimer = tic;
        coarseBaselineResult = fbtsSolve(coarseBaseCfg);
        fdtdTiming.sensitivity_baseline(transmitterIndex, iteration) = toc(fdtdTimer);
        coarsePerturbedCfg.source.samples(:) = 0;
        coarsePerturbedCfg.source.samples(transmitter, :) = coarseSourcePulse;
        phase = 'sensitivity_perturbation';
        fdtdTimer = tic;
        coarsePerturbedResult = fbtsSolve(coarsePerturbedCfg);
        fdtdTiming.sensitivity_perturbation(transmitterIndex, iteration) = toc(fdtdTimer);

        phaseTimer = tic;
        coarseSensitivityEz = ...
            (coarsePerturbedResult.rx_signals - ...
            coarseBaselineResult.rx_signals) ./ finiteDifferenceH;
        coarseSensitivityAll(transmitterIndex,:,:) = reshape( ...
            coarseSensitivityEz,1,numReceivers,coarseTemplate.Nt);
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

        phaseTiming.sensitivity_processing(transmitterIndex,iteration) = toc(phaseTimer);
        clear coarseBaselineResult coarsePerturbedResult
        fprintf('Sensitivity iteration %d/%d, transmitter %d/%d.\n', ...
            iteration, numIterations, transmitterIndex, numTransmitters);
    end

    phase = 'optimal_step_and_update'; phaseTimer = tic;
    %% Data-optimal step size, equations (31)-(33)
    assert(isfinite(stepA) && stepA > 0 && isfinite(stepQ), ...
        'fbts:InvalidOptimalStepSystem', ...
        'The optimal-step terms require finite q and finite positive a.');
    unconstrainedAlpha = stepQ ./ stepA;
    assert(isfinite(unconstrainedAlpha), 'fbts:InvalidOptimalStep', ...
        'The optimal permittivity step must be finite.');
    h.unconstrained_alpha(iteration) = unconstrainedAlpha;
    alpha = max(unconstrainedAlpha, 0);

    updatedEpsrDoi = s.epsrEst(doiMask) + ...
        alpha .* s.directionEpsr(doiMask);
    projectedCellCount = nnz(updatedEpsrDoi < epsrLowerBound);
    projectionApplied = projectedCellCount > 0;
    s.epsrEst(doiMask) = max(updatedEpsrDoi, epsrLowerBound);
    s.epsrEst(~doiMask) = epsrBackground(~doiMask);
    assert(all(isfinite(s.epsrEst),'all'),'fbts:NonfiniteEstimate','Updated material estimate is not finite.');

    s.previousProjectedGradEpsr = s.projectedGradEpsr;
    s.previousDirectionEpsr = s.directionEpsr;
    h.alpha(iteration) = alpha;
    h.pr_beta(iteration) = prBeta;
    h.pr_restart(iteration) = prRestart;
    h.projection_applied(iteration) = projectionApplied;
    h.projected_cell_count(iteration) = projectedCellCount;
    h.directional_derivative(iteration) = directionalDerivative;
    h.step_a(iteration) = stepA;
    h.step_q(iteration) = stepQ;
    h.finite_difference_h(iteration) = finiteDifferenceH;


    if execution.computeImageError
        h.relative_error_percent(iteration) = 100 .* ...
            norm(s.epsrEst(doiMask) - epsrTrue(doiMask)) ./ norm(epsrTrue(doiMask));
    else
        h.relative_error_percent(iteration) = NaN;
    end

    h.cost(iteration) = iterationCost;
    fprintf('Iteration %d cost: %.6e\n', iteration, iterationCost);
    fprintf(['Iteration %d optimal step: alpha = %.6e, ' ...
        'beta = %.6e, a = %.6e, q = %.6e, h = %.6e.\n'], ...
        iteration, alpha, prBeta, stepA, stepQ, finiteDifferenceH);
    if execution.computeImageError
        fprintf('Iteration %d/%d relative DOI error: %.6e%%.\n', ...
            iteration, numIterations, h.relative_error_percent(iteration));
    end
    h.iteration_runtime(iteration) = toc(iterationTimer);
    fprintf('Iteration %d/%d runtime: %.3f seconds.\n', ...
        iteration, numIterations, h.iteration_runtime(iteration));
    h.update_norm(iteration) = norm(s.epsrEst(doiMask)-previousEstimateDoi);
    s.lastModel = EzModel; s.lastSensitivity = EzSensitivity;
    phaseTiming.update(iteration) = toc(phaseTimer);
    phase = 'history_commit'; phaseTimer = tic;
    traceHistory.epsr_history_doi(:,iteration+1) = s.epsrEst(doiMask);
    if numIterations == 1
        % MAT files discard trailing singleton dimensions; matfile requires
        % exactly the number of stored dimensions when assigning a slice.
        traceHistory.rx_model_history(:,:,:) = EzModel;
        traceHistory.rx_sensitivity_coarse_history(:,:,:) = coarseSensitivityAll;
    else
        traceHistory.rx_model_history(:,:,:,iteration) = EzModel;
        traceHistory.rx_sensitivity_coarse_history(:,:,:,iteration) = coarseSensitivityAll;
    end
    traceHistory.complete_iterations = iteration;
    completedIterations = iteration;
    phaseTiming.history_write(iteration) = toc(phaseTimer);
    if ~isempty(execution.checkpointFile) && ...
            (mod(iteration,execution.checkpointEvery)==0 || iteration==numIterations)
        phase = 'checkpoint_commit';
        savedHistory = h;
        if ~execution.computeImageError, savedHistory = rmfield(savedHistory,'relative_error_percent'); end
        phaseTimer = tic;
        traceChecksum = historyChecksum(traceHistory,iteration,numIterations);
        phaseTiming.history_hash(iteration) = toc(phaseTimer);
        phaseTimer = tic;
        checkpoint = struct('identity',runIdentity,'completed_iterations',iteration, ...
            'state',s,'history',savedHistory,'fdtdTiming',fdtdTiming,'phaseTiming',phaseTiming,'trace_checksum',traceChecksum);
        checkpoint.checksum = fbtsHash(checkpoint);
        saveFbtsAtomic(execution.checkpointFile,struct('checkpoint',checkpoint));
        phaseTiming.checkpoint_write(iteration) = toc(phaseTimer);
    end
    phaseTiming.iteration_wall(iteration) = toc(iterationTimer);
    phase = 'iteration_callback';
    if ~isempty(execution.onIteration), execution.onIteration(iteration); end
end
catch exception
    if ~isempty(execution.checkpointFile)
        partialHistory = h;
        if ~execution.computeImageError, partialHistory = rmfield(partialHistory,'relative_error_percent'); end
        failureContext = struct('identity',runIdentity,'phase',phase,'iteration',iteration, ...
            'local_transmitter',transmitterIndex,'completed_iterations',completedIterations, ...
            'history',partialHistory,'timing',fdtdTiming,'current_estimate',s.epsrEst, ...
            'error_identifier',exception.identifier,'error_message',exception.message,'stack',exception.stack);
        try
            saveFbtsAtomic(fullfile(fileparts(execution.checkpointFile),'iteration_failure.mat'), ...
                struct('failureContext',failureContext));
        catch diagnostic
            warning('fbts:FailureDiagnosticSaveFailed','%s',diagnostic.message);
        end
    end
    rethrow(exception);
end

totalRuntime = sum(h.iteration_runtime);
averageRuntime = mean(h.iteration_runtime);
fprintf('Total FBTS runtime: %.3f seconds.\n', totalRuntime);
fprintf('Average FBTS iteration runtime: %.3f seconds.\n', averageRuntime);

%% Results
results = struct();
results.epsr_true = epsrTrue;
results.epsr_est = s.epsrEst;
results.gradient_epsr = s.projectedGradEpsr;
results.raw_gradient_epsr = s.gradEpsr;
results.direction_epsr = s.directionEpsr;
results.cost = h.cost;
results.alpha = h.alpha;
results.pr_beta = h.pr_beta;
results.pr_restart = h.pr_restart;
results.projection_applied = h.projection_applied;
results.projected_cell_count = h.projected_cell_count;
results.directional_derivative = h.directional_derivative;
results.epsr_lower_bound = epsrLowerBound;
results.step_a = h.step_a;
results.step_q = h.step_q;
results.finite_difference_h = h.finite_difference_h;
results.sensitivity_downsample_factor = sensitivityDownsampleFactor;
results.sensitivity_dt = h.sensitivity_dt;
results.sensitivity_Nt = h.sensitivity_Nt;
results.iteration_runtime = h.iteration_runtime;
results.total_runtime = totalRuntime;
results.average_runtime = averageRuntime;
results.iteration = (1:numIterations).';
results.relative_error_percent = h.relative_error_percent;
results.relative_error_definition = [ ...
    '100 * norm(epsr_est(doiMask) - epsr_true(doiMask)) / ' ...
    'norm(epsr_true(doiMask))'];
results.doi_mask = doiMask;
results.Ez_measured = EzMeasured;
results.Ez_model = s.lastModel;
results.Ez_sensitivity = s.lastSensitivity;

results.num_iterations = numIterations;
results.run_label = activeFbtsOptions.runLabel;
results.timing = struct('fdtd',fdtdTiming,'phases',phaseTiming);
results.history = rmfield(h,{'iteration_runtime','relative_error_percent'});
results.history.evaluated_state_index = (0:numIterations-1).';
results.history.updated_state_index = (1:numIterations).';
results.history.measurement_energy = sum(reshape(K,1,1,[]) .* EzMeasured.^2,'all')*cfg.dt;
results.history.channel_count = numTransmitters*numReceivers;
results.history.stopping_status = 'iteration_limit';
results.history.completed_iterations = numIterations;
results.epsr_history_doi = traceHistory.epsr_history_doi;
results.rx_model_history = traceHistory.rx_model_history;
results.rx_sensitivity_coarse_history = traceHistory.rx_sensitivity_coarse_history;
results.measurement_fingerprint = measurements.metadata.fingerprint;
results.code_hash = provenance.code_hash;
end

function checksum = historyChecksum(history,completed,requested)
if requested == 1
    model = history.rx_model_history(:,:,:);
    sensitivity = history.rx_sensitivity_coarse_history(:,:,:);
else
    model = history.rx_model_history(:,:,:,1:completed);
    sensitivity = history.rx_sensitivity_coarse_history(:,:,:,1:completed);
end
checksum = fbtsHash({history.epsr_history_doi(:,1:completed+1),model,sensitivity});
end

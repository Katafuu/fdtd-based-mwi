% run_dbim Run multi-frequency FDTD iterative Born imaging with a fixed analytical Green function.
% Run build_cfg_dbim.m first to construct the true synthetic cfg.
%
% Internal conventions:
%   material(x,y)
%   receiverFields(rx,tx,frequency)
%   doiFields(pixel,tx,frequency)

%% Repository / path setup
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end
workspaceRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(fullfile(workspaceRoot, 'buildLib'), '-end');
addpath(fullfile(workspaceRoot, 'forward_solver', 'mex'), '-end');
addpath(fullfile(scriptDir, 'lib'), '-begin');

assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'run_dbim:MissingConfig', ...
    'Run build_cfg_dbim.m before run_dbim.m.');
assert(exist('fdtd_mex', 'file') == 3, ...
    'run_dbim:MissingFdtdMex', ...
    'Build forward_solver/mex/fdtd_mex before running DBIM.');

%% Runtime options
defaultDbimOptions = struct( ...
    'frequencies', [0.9 1.0 1.1] .* 1e9, ...
    'maxIterations', 15, ...
    'relaxation', 0.5, ...
    'residualTolerance', 1e-3, ...
    'updateTolerance', 1e-3, ...
    'epsrBounds', [1 80], ...
    'sigmaBounds', [0 2], ...
    'epsrScale', 1.0, ...
    'sigmaScale', 0.1, ...
    'operatorPrecision', "single", ...
    'assemblyBlockSize', 64, ...
    'tsvdEnergy', 0.99, ...
    'tsvdInitialRank', 32, ...
    'tsvdMaximumRank', 512, ...
    'fullSvdMaxDimension', 512, ...
    'outputDirectory', fullfile(scriptDir, 'figs'), ...
    'plotResolution', 300, ...
    'enablePlot', true, ...
    'verbose', true);
if ~exist('dbimOpts', 'var') || isempty(dbimOpts)
    dbimOpts = struct();
end
dbimOpts = mergeDbimOptions(dbimOpts, defaultDbimOptions);
validateDbimOptions(dbimOpts);
frequencies = double(dbimOpts.frequencies(:).');
runOutputDirectory = createNextRunDirectory(dbimOpts.outputDirectory);
totalRuntimeTimer = tic;

%% True medium, background estimate, DOI, and reciprocal data pairs
trueCfg = cfg;
cfg_est = fdtdmat.setBackgroundDefault(trueCfg, trueCfg.grid.background);
if isfield(cfg_est, 'targets')
    cfg_est.targets = struct([]);
end

doiMask = logical(trueCfg.antennas.doiMask);
doiLinearIndex = find(doiMask);
[doiX, doiY] = ind2sub([trueCfg.Nx trueCfg.Ny], doiLinearIndex);
doiPositions = [(doiX - 1) .* trueCfg.dx, ...
    (doiY - 1) .* trueCfg.dy];
antennaPositions = (double(trueCfg.antennas.pos) - 1) .* ...
    [trueCfg.dx trueCfg.dy];

numAntennas = trueCfg.antennas.numAntennas;
[pairTx, pairRx] = find(triu(true(numAntennas), 1));
pairs = [pairTx pairRx];
numDoiPixels = numel(doiLinearIndex);

%% Fixed analytical homogeneous Green function
background = cfg_est.grid.background;
numFrequencies = numel(frequencies);
greenWavenumber = complex(zeros(1, numFrequencies));
analyticGreen = complex(zeros(numDoiPixels, numAntennas, ...
    numFrequencies, char(dbimOpts.operatorPrecision)));
for frequencyIndex = 1:numFrequencies
    omega = 2*pi*frequencies(frequencyIndex);
    relativeComplexPermittivity = background.epsr - ...
        1i*background.cond_e/(omega*cfg_est.eps0);
    greenWavenumber(frequencyIndex) = omega*sqrt( ...
        cfg_est.mu0*cfg_est.eps0*relativeComplexPermittivity);
    analyticGreen(:, :, frequencyIndex) = green2D( ...
        doiPositions, antennaPositions, greenWavenumber(frequencyIndex));
end
greenFunctionInfo = struct( ...
    'type', "analytical homogeneous 2-D", ...
    'wavenumber', greenWavenumber, ...
    'backgroundEpsr', background.epsr, ...
    'backgroundSigma', background.cond_e);

fprintf('FDTD-DBIM: %d DOI pixels, %d antenna pairs, %d frequencies.\n', ...
    numDoiPixels, size(pairs, 1), numel(frequencies));

%% Fixed synthetic measurement data
fprintf('Generating fixed synthetic target measurements.\n');
measurementTimer = tic;
measurementScan = runFdtdDbimScan( ...
    trueCfg, frequencies, doiLinearIndex, false);
measurementRuntime = toc(measurementTimer);
measuredFields = measurementScan.receiverFields;
sourceSpectrum = measurementScan.sourceSpectrum;
clear measurementScan
fprintf('Synthetic measurement runtime: %.3f s.\n', measurementRuntime);

%% DBIM iteration state
maxEvaluations = dbimOpts.maxIterations + 1;
evaluatedIteration = nan(maxEvaluations, 1);
residualHistory = nan(maxEvaluations, 1);
updateNormHistory = nan(maxEvaluations, 1);
tsvdRankHistory = nan(maxEvaluations, 1);
tsvdEnergyHistory = nan(maxEvaluations, 1);
forwardTimeHistory = nan(maxEvaluations, 1);
solveTimeHistory = nan(maxEvaluations, 1);
numTimingVariables = 4;
computationTime = nan(numTimingVariables, dbimOpts.maxIterations);
iterationRuntime = nan(1, dbimOpts.maxIterations);
computationTimeLabels = [
    "Background/current-medium FDTD simulation"
    "FFT / frequency extraction"
    "Assembly of A and b"
    "Linear-system solve"
];
iterationsRun = 0;

evaluationCount = 0;
updatesApplied = 0;
bestResidual = Inf;
bestIteration = 0;
bestCfg = cfg_est;
bestPredictedFields = complex(zeros(size(measuredFields)));
stopReason = "maximum iterations";

linearSystemOptions = struct( ...
    'epsrScale', dbimOpts.epsrScale, ...
    'sigmaScale', dbimOpts.sigmaScale, ...
    'operatorPrecision', dbimOpts.operatorPrecision, ...
    'assemblyBlockSize', dbimOpts.assemblyBlockSize);
tsvdOptions = struct( ...
    'energyThreshold', dbimOpts.tsvdEnergy, ...
    'initialRank', dbimOpts.tsvdInitialRank, ...
    'maximumRank', dbimOpts.tsvdMaximumRank, ...
    'fullSvdMaxDimension', dbimOpts.fullSvdMaxDimension);
updateOptions = struct( ...
    'relaxation', dbimOpts.relaxation, ...
    'epsrBounds', dbimOpts.epsrBounds, ...
    'sigmaBounds', dbimOpts.sigmaBounds, ...
    'epsrScale', dbimOpts.epsrScale, ...
    'sigmaScale', dbimOpts.sigmaScale);

%% Born outer loop with a fixed analytical Green function
for iteration = 1:dbimOpts.maxIterations
    iterationTimer = tic;
    if dbimOpts.verbose
        fprintf('DBIM iteration %d/%d: current-medium FDTD scan.\n', ...
            iteration, dbimOpts.maxIterations);
    end
    currentScan = runFdtdDbimScan( ...
        cfg_est, frequencies, doiLinearIndex, true);
    iterationsRun = iteration;
    computationTime(1, iteration) = currentScan.simulationRuntime;
    computationTime(2, iteration) = currentScan.fourierRuntime;


    assemblyTimer = tic;
    [operator, rhs, systemMetadata] = buildDbimLinearSystem( ...
        measuredFields, currentScan.receiverFields, currentScan.doiFields, ...
        pairs, analyticGreen, frequencies, cfg_est, ...
        linearSystemOptions);
    computationTime(3, iteration) = toc(assemblyTimer);

    evaluationCount = evaluationCount + 1;
    evaluatedIteration(evaluationCount) = updatesApplied;
    residualHistory(evaluationCount) = systemMetadata.normalizedResidual;
    forwardTimeHistory(evaluationCount) = currentScan.runtime;

    if systemMetadata.normalizedResidual < bestResidual
        bestResidual = systemMetadata.normalizedResidual;
        bestIteration = updatesApplied;
        bestCfg = cfg_est;
        bestPredictedFields = currentScan.receiverFields;
    end

    fprintf('DBIM iteration %d: normalized residual %.6e.\n', ...
        iteration, systemMetadata.normalizedResidual);
    if systemMetadata.normalizedResidual <= dbimOpts.residualTolerance
        stopReason = "residual tolerance";
        iterationRuntime(iteration) = toc(iterationTimer);
        fprintf('DBIM iteration %d total runtime: %.3f s.\n', ...
            iteration, iterationRuntime(iteration));
        clear operator rhs currentScan
        break
    end

    solveTimer = tic;
    [scaledUpdate, tsvdInfo] = solve_tsvd(operator, rhs, tsvdOptions);
    computationTime(4, iteration) = toc(solveTimer);
    solveTimeHistory(evaluationCount) = computationTime(4, iteration);
    tsvdRankHistory(evaluationCount) = tsvdInfo.retainedRank;
    tsvdEnergyHistory(evaluationCount) = tsvdInfo.achievedEnergy;
    if tsvdInfo.hitRankCap
        warning('run_dbim:TsvdRankCap', ...
            'TSVD reached rank %d with %.3f%% of operator energy.', ...
            tsvdInfo.retainedRank, 100*tsvdInfo.achievedEnergy);
    end

    deltaEpsr = dbimOpts.epsrScale .* scaledUpdate(1:numDoiPixels);
    deltaSigma = dbimOpts.sigmaScale .* ...
        scaledUpdate(numDoiPixels + (1:numDoiPixels));
    [cfg_est, updateInfo] = applyDbimUpdate( ...
        cfg_est, deltaEpsr, deltaSigma, doiMask, updateOptions);
    updatesApplied = updatesApplied + 1;
    updateNormHistory(evaluationCount) = updateInfo.relativeUpdateNorm;

    fprintf(['DBIM iteration %d: rank %d, energy %.3f%%, ' ...
        'relative update %.6e.\n'], iteration, tsvdInfo.retainedRank, ...
        100*tsvdInfo.achievedEnergy, updateInfo.relativeUpdateNorm);
    if dbimOpts.verbose
        fprintf(['DBIM iteration %d timing [s]: FDTD %.3f, FFT %.3f, ' ...
            'assembly %.3f, solve %.3f.\n'], iteration, ...
            computationTime(:, iteration));
    end
    iterationRuntime(iteration) = toc(iterationTimer);
    fprintf('DBIM iteration %d total runtime: %.3f s.\n', ...
        iteration, iterationRuntime(iteration));

    clear operator rhs currentScan scaledUpdate deltaEpsr deltaSigma
    if updateInfo.relativeUpdateNorm <= dbimOpts.updateTolerance
        stopReason = "update tolerance";
        break
    end
end

%% Evaluate the last applied update and return the best-residual estimate
fprintf('Evaluating the final DBIM material estimate.\n');
finalCandidateScan = runFdtdDbimScan( ...
    cfg_est, frequencies, doiLinearIndex, false);
[finalCandidateResidual, ~] = normalizedReceiverResidual( ...
    measuredFields, finalCandidateScan.receiverFields, pairs);

evaluationCount = evaluationCount + 1;
evaluatedIteration(evaluationCount) = updatesApplied;
residualHistory(evaluationCount) = finalCandidateResidual;
forwardTimeHistory(evaluationCount) = finalCandidateScan.runtime;

if true
    bestResidual = finalCandidateResidual;
    bestIteration = updatesApplied;
    bestCfg = cfg_est;
    bestPredictedFields = finalCandidateScan.receiverFields;
end

cfg_est = bestCfg;
finalPredictedFields = bestPredictedFields;
totalRuntime = toc(totalRuntimeTimer);

%% Workspace result and plots
historyIndex = 1:evaluationCount;
computationTime = computationTime(:, 1:iterationsRun);
iterationRuntime = iterationRuntime(1:iterationsRun);
dbimResult = struct();
dbimResult.estimate = struct( ...
    'epsr', cfg_est.grid.epsr, ...
    'sigma', cfg_est.grid.cond_e);
dbimResult.truth = struct( ...
    'epsr', trueCfg.grid.epsr, ...
    'sigma', trueCfg.grid.cond_e);
dbimResult.options = dbimOpts;
dbimResult.frequencies = frequencies;
dbimResult.pairs = pairs;
dbimResult.doiMask = doiMask;
dbimResult.greenFunction = greenFunctionInfo;
dbimResult.sourceSpectrum = sourceSpectrum;
dbimResult.measuredFields = measuredFields;
dbimResult.predictedFields = finalPredictedFields;
dbimResult.finalResidual = bestResidual;
dbimResult.finalCandidateResidual = finalCandidateResidual;
dbimResult.bestIteration = bestIteration;
dbimResult.iterationsCompleted = updatesApplied;
dbimResult.stopReason = stopReason;
dbimResult.computationTime = computationTime;
dbimResult.computationTimeLabels = computationTimeLabels;
dbimResult.measurementRuntime = measurementRuntime;
dbimResult.iterationRuntime = iterationRuntime;
dbimResult.totalRuntime = totalRuntime;
dbimResult.history = struct( ...
    'evaluatedIteration', evaluatedIteration(historyIndex), ...
    'residual', residualHistory(historyIndex), ...
    'updateNorm', updateNormHistory(historyIndex), ...
    'tsvdRank', tsvdRankHistory(historyIndex), ...
    'tsvdEnergy', tsvdEnergyHistory(historyIndex), ...
    'forwardTime', forwardTimeHistory(historyIndex), ...
    'solveTime', solveTimeHistory(historyIndex));

fprintf(['FDTD-DBIM complete: %d updates, best iteration %d, ' ...
    'normalized residual %.6e (%s).\n'], updatesApplied, ...
    bestIteration, bestResidual, stopReason);

figureFiles = plotDbimResults( ...
    trueCfg, dbimResult, runOutputDirectory, ...
    dbimOpts.enablePlot, dbimOpts.plotResolution);
resultFile = string(fullfile( ...
    runOutputDirectory, 'dbim_run_data.mat'));
dbimResult.outputFiles = struct( ...
    'figures', figureFiles, ...
    'matFile', resultFile);
save(char(resultFile), 'cfg', 'cfg_est', 'trueCfg', 'dbimOpts', ...
    'computationTime', 'computationTimeLabels', 'measurementRuntime', ...
    'iterationRuntime', 'totalRuntime', 'dbimResult', '-v7.3');
fprintf('Saved DBIM figures and run data under %s.\n', ...
    runOutputDirectory);
fprintf('Total DBIM computation runtime: %.3f s.\n', totalRuntime);

function runDirectory = createNextRunDirectory(outputRoot)
outputRoot = char(outputRoot);
if ~isfolder(outputRoot)
    [created, message] = mkdir(outputRoot);
    if ~created
        error('run_dbim:CreateOutputRootFailed', '%s', message);
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
    error('run_dbim:CreateRunDirectoryFailed', '%s', message);
end
end

function options = mergeDbimOptions(options, defaults)
if ~isstruct(options) || ~isscalar(options)
    error('run_dbim:InvalidOptions', 'dbimOpts must be a scalar struct.');
end
defaultNames = fieldnames(defaults);
for optionIndex = 1:numel(defaultNames)
    optionName = defaultNames{optionIndex};
    if ~isfield(options, optionName) || isempty(options.(optionName))
        options.(optionName) = defaults.(optionName);
    end
end
end

function validateDbimOptions(options)
if ~isnumeric(options.frequencies) || isempty(options.frequencies) || ...
        any(~isfinite(options.frequencies)) || any(options.frequencies <= 0)
    error('run_dbim:InvalidFrequencies', ...
        'dbimOpts.frequencies must contain finite positive values.');
end
positiveIntegerNames = {'maxIterations', 'assemblyBlockSize', ...
    'tsvdInitialRank', 'tsvdMaximumRank', 'fullSvdMaxDimension', ...
    'plotResolution'};
for optionIndex = 1:numel(positiveIntegerNames)
    value = options.(positiveIntegerNames{optionIndex});
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
            value < 1 || value ~= round(value)
        error('run_dbim:InvalidOptions', ...
            'dbimOpts.%s must be a positive integer.', ...
            positiveIntegerNames{optionIndex});
    end
end
unitIntervalNames = {'relaxation', 'tsvdEnergy'};
for optionIndex = 1:numel(unitIntervalNames)
    value = options.(unitIntervalNames{optionIndex});
    if ~isnumeric(value) || ~isscalar(value) || ...
            ~isfinite(value) || value <= 0 || value > 1
        error('run_dbim:InvalidOptions', ...
            'dbimOpts.%s must lie in (0,1].', unitIntervalNames{optionIndex});
    end
end
if ~isscalar(options.enablePlot) || ~isscalar(options.verbose)
    error('run_dbim:InvalidOptions', ...
        'dbimOpts.enablePlot and verbose must be scalars.');
end
if ~(ischar(options.outputDirectory) || ...
        (isstring(options.outputDirectory) && ...
        isscalar(options.outputDirectory))) || ...
        strlength(string(options.outputDirectory)) == 0
    error('run_dbim:InvalidOptions', ...
        'dbimOpts.outputDirectory must be a nonempty text scalar.');
end
end

function [normalizedResidual, complexResidual] = normalizedReceiverResidual( ...
    measuredFields, predictedFields, pairs)
numAntennas = size(measuredFields, 1);
numFrequencies = size(measuredFields, 3);
pairIndex = sub2ind([numAntennas numAntennas], pairs(:, 2), pairs(:, 1));
complexResidual = complex(zeros(size(pairs, 1), numFrequencies));
measuredPairs = complex(zeros(size(pairs, 1), numFrequencies));
for frequencyIndex = 1:numFrequencies
    measuredAtFrequency = measuredFields(:, :, frequencyIndex);
    predictedAtFrequency = predictedFields(:, :, frequencyIndex);
    measuredPairs(:, frequencyIndex) = measuredAtFrequency(pairIndex);
    complexResidual(:, frequencyIndex) = ...
        measuredAtFrequency(pairIndex) - predictedAtFrequency(pairIndex);
end
normalizedResidual = norm(complexResidual(:)) / ...
    max(norm(measuredPairs(:)), eps);
end

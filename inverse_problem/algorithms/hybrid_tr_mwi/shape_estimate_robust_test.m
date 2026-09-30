% Run three independent robustness cases using the shape_estimate workflow.
% Each case shares one target geometry and saves figures plus MATLAB data.

%% Editable experiment controls
SNR_dB = 20;
targetEpsrBase = 4.0;
targetEpsrVariance = 0.04;
backgroundEpsrBase = 1.05;
backgroundEpsrVariance = 0.0004;
inhomCorrelationCells = 6;
robustSeed = 31415;

validateattributes(SNR_dB, {'numeric'}, {'real', 'scalar', 'finite'});
validateattributes(targetEpsrBase, {'numeric'}, ...
    {'real', 'scalar', 'finite', '>', 1});
validateattributes(backgroundEpsrBase, {'numeric'}, ...
    {'real', 'scalar', 'finite', '>=', 1});
validateattributes(targetEpsrVariance, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'nonnegative'});
validateattributes(backgroundEpsrVariance, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'nonnegative'});
validateattributes(inhomCorrelationCells, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'positive'});
validateattributes(robustSeed, {'numeric'}, ...
    {'real', 'scalar', 'integer', '>=', 0, '<=', 2^32-4});

%% Build one shared geometry and nominal material configuration
scriptDir = fileparts(mfilename('fullpath'));
rng(robustSeed, 'twister');
run(fullfile(scriptDir, 'build_cfg.m'));
assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'shape_estimate_robust_test:MissingConfig', ...
    'build_cfg.m did not create a configuration.');
assert(exist('fdtd_mex', 'file') == 3, ...
    'shape_estimate_robust_test:MissingFdtdMex', ...
    'Build forward_solver/mex/fdtd_mex before running this script.');

cfg.opts.numIterations = 1;
cfg.opts.storeFieldHistory = false; % The saved focus matrices do not need Ez history.

% Capture the target before changing either material map. Background
% variation must not turn the entire varied region into a target outline.
epsrTolerance = 10*eps(max(1, max(abs(cfg.grid.epsr(:)))));
conductivityTolerance = 10*eps(max(1, max(abs(cfg.grid.cond_e(:)))));
targetMask = ...
    abs(cfg.grid.epsr - cfg.grid.background.epsr) > epsrTolerance | ...
    abs(cfg.grid.cond_e - cfg.grid.background.cond_e) > conductivityTolerance;
assert(any(targetMask(:)), ...
    'shape_estimate_robust_test:EmptyTarget', ...
    'The shared configuration contains no target cells.');

% Set the same nominal target and background in all three cases.
cfg.grid.background.epsr = backgroundEpsrBase;
cfg.grid.epsr(:) = backgroundEpsrBase;
cfg.grid.epsr(targetMask) = targetEpsrBase;
cfg.grid.epsr_bg(:) = backgroundEpsrBase;
for targetIdx = 1:numel(cfg.targets)
    cfg.targets(targetIdx).material.epsr = targetEpsrBase;
end

% The PML conductivity depends on background impedance. Rebuild its profile
% after changing the nominal background permittivity.
cfg.pml.eta0 = sqrt((cfg.mu0 * cfg.grid.background.murx) / ...
    (cfg.eps0 * backgroundEpsrBase));
cfg.pml.sigma_max = -(cfg.pml.m + 1) * log(cfg.pml.R) / ...
    (2 * cfg.pml.eta0 * cfg.pml.physicalThickness);
pmlProfile = fdtdpml.build_rectangularPML(cfg.grid, ...
    cfg.pml.thickness, cfg.pml.m, cfg.pml.sigma_max, cfg.pml.kappaMax);
profileFields = fieldnames(pmlProfile);
for fieldIdx = 1:numel(profileFields)
    fieldName = profileFields{fieldIdx};
    cfg.pml.(fieldName) = pmlProfile.(fieldName);
end
cfg.pml.sigma_e = max(cfg.pml.condx, cfg.pml.condy);

baselineCfg = cfg;
backgroundMask = ~targetMask & ~cfg.pml.mask;

%% Find a new run folder without overwriting previous results
figsRoot = fullfile(scriptDir, 'figs', 'shape_estimate');
if ~isfolder(figsRoot)
    mkdir(figsRoot);
end
existingRuns = dir(fullfile(figsRoot, 'run_*'));
existingRunNumbers = zeros(0, 1);
for directoryIdx = 1:numel(existingRuns)
    if ~existingRuns(directoryIdx).isdir
        continue;
    end
    token = regexp(existingRuns(directoryIdx).name, ...
        '^run_(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        existingRunNumbers(end + 1, 1) = str2double(token{1}); %#ok<SAGROW>
    end
end
if isempty(existingRunNumbers)
    nextRunNumber = 1;
else
    nextRunNumber = max(existingRunNumbers) + 1;
end
runDir = fullfile(figsRoot, sprintf('run_%04d', nextRunNumber));
mkdir(runDir);

%% Run the three cases separately
caseNames = {'noise', 'inhom_targ', 'inhom_bg'};
for caseIdx = 1:numel(caseNames)
    caseName = caseNames{caseIdx};
    cfg = baselineCfg;
    cfg.opts.SNR_dB = [];
    inhomParams = struct();
    inhomParams.target = struct( ...
        'enabled', false, 'baseEpsr', targetEpsrBase, ...
        'variance', targetEpsrVariance, ...
        'correlationCells', inhomCorrelationCells, ...
        'seed', robustSeed + 2, 'realizedMean', NaN, ...
        'realizedVariance', NaN, 'minEpsr', NaN, 'maxEpsr', NaN);
    inhomParams.background = struct( ...
        'enabled', false, 'baseEpsr', backgroundEpsrBase, ...
        'variance', backgroundEpsrVariance, ...
        'correlationCells', inhomCorrelationCells, ...
        'seed', robustSeed + 3, 'realizedMean', NaN, ...
        'realizedVariance', NaN, 'minEpsr', NaN, 'maxEpsr', NaN);

    switch caseName
        case 'noise'
            cfg.opts.SNR_dB = SNR_dB;
            rng(robustSeed + 1, 'twister');
        case 'inhom_targ'
            rng(inhomParams.target.seed, 'twister');
            targetMap = makeSmoothEpsrMap( ...
                targetMask, targetEpsrBase, targetEpsrVariance, ...
                inhomCorrelationCells);
            cfg.grid.epsr(targetMask) = targetMap(targetMask);
            inhomParams.target = withRealizedStats( ...
                inhomParams.target, cfg.grid.epsr(targetMask));
        case 'inhom_bg'
            rng(inhomParams.background.seed, 'twister');
            backgroundMap = makeSmoothEpsrMap( ...
                backgroundMask, backgroundEpsrBase, ...
                backgroundEpsrVariance, inhomCorrelationCells);
            cfg.grid.epsr(backgroundMask) = backgroundMap(backgroundMask);
            inhomParams.background = withRealizedStats( ...
                inhomParams.background, cfg.grid.epsr(backgroundMask));
    end

    figsDir = fullfile(runDir, caseName);
    mkdir(figsDir);
    focusPoints = zeros(cfg.antennas.numAntennas, 2);
    focusMagFrames = zeros(cfg.Nx, cfg.Ny, cfg.antennas.numAntennas);
    focusMagnitudeImages = zeros(size(focusMagFrames));
    receiverNoisePowers = nan(cfg.antennas.numAntennas, 1);
    totalRuntime = 0;

    fprintf('\nRunning %s case under %s.\n', caseName, figsDir);
    for antennaIdx = 1:cfg.antennas.numAntennas
        runTimer = tic;
        antennaCfg = cfg;
        antennaCfg.antennas.txAntennas = antennaIdx;
        antennaCfg.source.samples(:) = 0;
        antennaCfg.source.samples(antennaIdx, :) = ...
            antennaCfg.source.func((0:antennaCfg.Nt-1) .* antennaCfg.dt);

        if strcmp(caseName, 'noise')
            itrResult = itr_run(antennaCfg, SNR_dB);
            [~, focusLinearIndex] = max(abs(itrResult.focusMagFrame(:)));
            [focusX, focusY] = ind2sub( ...
                size(itrResult.focusMagFrame), focusLinearIndex);
            focusPoint = [focusX, focusY];
            receiverNoisePowers(antennaIdx) = ...
                itrResult.receiverNoisePowers(end);
        else
            [itrResult, focusPoint] = tr_mwi(antennaCfg);
        end

        focusPoints(antennaIdx, :) = focusPoint;
        focusMagFrames(:, :, antennaIdx) = itrResult.focusMagFrame;
        focusMagnitudeImages(:, :, antennaIdx) = ...
            abs(itrResult.focusMagFrame);
        saveMaxEz(itrResult, focusPoint, antennaIdx, ...
            figsDir, targetMask, cfg);
        clear itrResult;

        elapsedTime = toc(runTimer);
        totalRuntime = totalRuntime + elapsedTime;
        fprintf('Case %s: antenna %d/%d done in %.2f seconds.\n', ...
            caseName, antennaIdx, cfg.antennas.numAntennas, elapsedTime);
    end

    fprintf('Case %s total runtime: %.2f seconds.\n', caseName, totalRuntime);
    fprintf('Average runtime per antenna: %.2f seconds.\n', ...
        totalRuntime/cfg.antennas.numAntennas);
    saveFocusOverview(focusPoints, cfg, targetMask, figsDir);

    % Save the numeric fields represented by the PNGs, shape points, the
    % complete case configuration, and nominal/realized material settings.
    [doiPoints, doiPointIndices] = selectShapePointsInsideDoi( ...
        focusPoints, cfg.antennas.doiMask);
    [inlierPoints, ~, inlierDoiIndices] = ...
        filterShapeOutlinePoints(doiPoints, 1.5);
    estimatedOutlinePoints = orderShapeOutlinePoints( ...
        inlierPoints, doiPointIndices(inlierDoiIndices));
    antennaIndices = (1:cfg.antennas.numAntennas).';
    caseSNR_dB = cfg.opts.SNR_dB;
    save(fullfile(figsDir, 'shape_estimate_results.mat'), ...
        'cfg', 'caseName', 'caseSNR_dB', 'robustSeed', ...
        'inhomParams', 'targetMask', 'backgroundMask', ...
        'antennaIndices', 'focusPoints', 'estimatedOutlinePoints', ...
        'focusMagFrames', 'focusMagnitudeImages', ...
        'receiverNoisePowers', '-v7.3');
end

fprintf('\nSaved all robustness cases under %s.\n', runDir);

%% Plotting functions copied from shape_estimate.m
function saveMaxEz(itrResult, focusPoint, antennaIdx, figsDir, targetMask, cfg)
    maxMagnitudeFigure = figure('Visible', 'off');
    axesHandle = axes(maxMagnitudeFigure);
    imagesc(axesHandle, abs(itrResult.focusMagFrame).');
    axis(axesHandle, 'equal', 'tight');
    set(axesHandle, 'YDir', 'normal');
    colorbar(axesHandle);
    hold(axesHandle, 'on');
    contour(axesHandle, double(cfg.antennas.doiMask).', [0.5 0.5], ...
        'Color', [0 0.9 0.9], 'LineStyle', '--', 'LineWidth', 1.5);
    contour(axesHandle, double(targetMask).', [0.5 0.5], ...
        'Color', [1 0.25 0.25], 'LineWidth', 1.75);
    plot(axesHandle, cfg.antennas.pos(:, 1), cfg.antennas.pos(:, 2), ...
        'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
    plot(axesHandle, cfg.antennas.pos(antennaIdx, 1), ...
        cfg.antennas.pos(antennaIdx, 2), 'rx', ...
        'MarkerSize', 9, 'LineWidth', 1.5);
    plot(axesHandle, focusPoint(1), focusPoint(2), 'wo', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.5);
    hold(axesHandle, 'off');
    xlabel(axesHandle, 'x grid index');
    ylabel(axesHandle, 'y grid index');
    title(axesHandle, ...
        sprintf('Maximum-magnitude reconstruction, antenna %d', antennaIdx));
    exportgraphics(maxMagnitudeFigure, ...
        fullfile(figsDir, sprintf('%d.png', antennaIdx)), ...
        'Resolution', 300);
    close(maxMagnitudeFigure);
end

function saveFocusOverview(points, cfg, targetMask, figsDir)
    [pointsFigure, outlineFigure] = ...
        plotShapeEstimateOverview(points, cfg, targetMask, 1.5);
    figureCleanup = onCleanup( ...
        @() close([pointsFigure outlineFigure]));

    exportgraphics(pointsFigure, ...
        fullfile(figsDir, 'all_maximum_points.png'), 'Resolution', 300);
    exportgraphics(outlineFigure, ...
        fullfile(figsDir, 'estimated_outline.png'), 'Resolution', 300);
end

function epsrMap = makeSmoothEpsrMap(mask, baseEpsr, variance, correlationCells)
    epsrMap = baseEpsr * ones(size(mask));
    if variance == 0
        return;
    end

    kernelRadius = ceil(3 * correlationCells);
    kernelAxis = -kernelRadius:kernelRadius;
    kernel = exp(-0.5 * (kernelAxis ./ correlationCells).^2);
    kernel = kernel ./ sum(kernel);
    smoothField = conv2(kernel(:), kernel(:).', ...
        randn(size(mask)), 'same');
    maskedValues = smoothField(mask);
    maskedStd = std(maskedValues);
    if maskedStd == 0
        error('shape_estimate_robust_test:ConstantRandomField', ...
            'The random field has no variation within the selected region.');
    end
    normalizedField = (smoothField - mean(maskedValues)) ./ maskedStd;
    normalizedField = max(-2.5, min(2.5, normalizedField));
    epsrMap = max(1.0, baseEpsr + sqrt(variance) .* normalizedField);
end

function params = withRealizedStats(params, values)
    params.enabled = true;
    params.realizedMean = mean(values);
    params.realizedVariance = var(values);
    params.minEpsr = min(values);
    params.maxEpsr = max(values);
end

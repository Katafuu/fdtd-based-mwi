% Run three independent robustness cases for each target setup.
% Each setup gets its own run folder with noise, target, and background cases.

%% Editable experiment controls
if ~exist('setups', 'var')
    setups = {struct()}; % Add struct entries to run more target setups.
end
SNR_dB = 20;
targetEpsrVariance = 0.04;
backgroundEpsrBase = 1.05;
backgroundEpsrVariance = 0.0004;
inhomCorrelationCells = 6;
robustSeed = 31415;

assert(iscell(setups) && ~isempty(setups) && ...
    all(cellfun(@(item) isstruct(item) && isscalar(item), setups)), ...
    'shape_estimate_robust_test:InvalidSetups', ...
    'setups must be a nonempty cell array of scalar setup structs.');
validateattributes(SNR_dB, {'numeric'}, {'real', 'scalar', 'finite'});
validateattributes(backgroundEpsrBase, {'numeric'}, ...
    {'real', 'scalar', 'finite', '>=', 1});
validateattributes(targetEpsrVariance, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'nonnegative'});
validateattributes(backgroundEpsrVariance, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'nonnegative'});
validateattributes(inhomCorrelationCells, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'positive'});
validateattributes(robustSeed, {'numeric'}, ...
    {'real', 'scalar', 'integer', '>=', 0, '<=', 2^32-1});

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir, '-begin');
assert(exist('fdtd_mex', 'file') == 3 || ...
    isfile(fullfile(scriptDir, '..', '..', '..', 'forward_solver', ...
    'mex', ['fdtd_mex.' mexext])), ...
    'shape_estimate_robust_test:MissingFdtdMex', ...
    'Build forward_solver/mex/fdtd_mex before running this script.');

%% Find the next available run number
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

caseNames = {'noise', 'inhom_targ', 'inhom_bg'};
for setupIdx = 1:numel(setups)
    setup = setups{setupIdx};
    cfg = build_cfg(setup);
    cfg.opts.numIterations = 1;
    cfg.opts.storeFieldHistory = false;

    % Include the coarse target masks and any boundary cells retained by
    % the high-resolution material rasterization.
    targetMask = false(cfg.Nx, cfg.Ny);
    for targetIdx = 1:numel(cfg.targets)
        targetMask = targetMask | cfg.targets(targetIdx).mask;
    end
    epsrTolerance = 10 * eps(max(1, max(abs(cfg.grid.epsr(:)))));
    condTolerance = 10 * eps(max(1, max(abs(cfg.grid.cond_e(:)))));
    targetMask = targetMask | ...
        abs(cfg.grid.epsr - cfg.grid.background.epsr) > epsrTolerance | ...
        abs(cfg.grid.cond_e - cfg.grid.background.cond_e) > condTolerance;
    assert(any(targetMask(:)), ...
        'shape_estimate_robust_test:EmptyTarget', ...
        'The setup contains no target cells.');

    % Give all three cases the same nominal background and update the PML
    % impedance to match it. Existing target permittivities are preserved.
    cfg.grid.background.epsr = backgroundEpsrBase;
    cfg.grid.epsr(~targetMask) = backgroundEpsrBase;
    cfg.grid.epsr_bg(:) = backgroundEpsrBase;
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

    runDir = fullfile(figsRoot, sprintf('run_%04d', nextRunNumber));
    nextRunNumber = nextRunNumber + 1;
    mkdir(runDir);

    %% Run the three cases separately
    for caseIdx = 1:numel(caseNames)
        caseName = caseNames{caseIdx};
        caseSeed = mod(double(robustSeed) + ...
            3 * (setupIdx - 1) + caseIdx, 2^32);
        cfg = baselineCfg;
        cfg.opts.SNR_dB = [];
        info = struct();
        switch caseName
            case 'noise'
                cfg.opts.SNR_dB = SNR_dB;
                rng(caseSeed, 'twister');
            case 'inhom_targ'
                [cfg, info] = fdtdmat.applyEpsrInhomogeneity( ...
                    cfg, targetMask, struct( ...
                    'variance', targetEpsrVariance, ...
                    'correlationCells', inhomCorrelationCells, ...
                    'seed', caseSeed));
            case 'inhom_bg'
                [cfg, info] = fdtdmat.applyEpsrInhomogeneity( ...
                    cfg, backgroundMask, struct( ...
                    'variance', backgroundEpsrVariance, ...
                    'correlationCells', inhomCorrelationCells, ...
                    'seed', caseSeed));
        end

        figsDir = fullfile(runDir, caseName);
        mkdir(figsDir);
        focusPoints = zeros(cfg.antennas.numAntennas, 2);
        focusMagFrames = zeros(cfg.Nx, cfg.Ny, cfg.antennas.numAntennas);
        focusMagnitudeImages = zeros(size(focusMagFrames));
        receiverNoisePowers = nan(cfg.antennas.numAntennas, 1);
        totalRuntime = 0;

        fprintf('\nSetup %d/%d: running %s under %s.\n', ...
            setupIdx, numel(setups), caseName, figsDir);
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

        fprintf('Case %s total runtime: %.2f seconds.\n', ...
            caseName, totalRuntime);
        fprintf('Average runtime per antenna: %.2f seconds.\n', ...
            totalRuntime/cfg.antennas.numAntennas);
        saveFocusOverview(focusPoints, cfg, targetMask, figsDir);

        [doiPoints, doiPointIndices] = selectShapePointsInsideDoi( ...
            focusPoints, cfg.antennas.doiMask);
        [inlierPoints, ~, inlierDoiIndices] = ...
            filterShapeOutlinePoints(doiPoints, 1.5);
        estimatedOutlinePoints = orderShapeOutlinePoints( ...
            inlierPoints, doiPointIndices(inlierDoiIndices));
        antennaIndices = (1:cfg.antennas.numAntennas).';
        caseSNR_dB = cfg.opts.SNR_dB;
        save(fullfile(figsDir, 'shape_estimate_results.mat'), ...
            'cfg', 'setup', 'setupIdx', 'caseName', 'caseSeed', ...
            'caseSNR_dB', 'robustSeed', 'info', ...
            'targetMask', 'backgroundMask', ...
            'antennaIndices', 'focusPoints', 'estimatedOutlinePoints', ...
            'focusMagFrames', 'focusMagnitudeImages', ...
            'receiverNoisePowers', '-v7.3');
    end
    fprintf('\nSaved setup %d under %s.\n', setupIdx, runDir);
end

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

% Run one time-reversal iteration per antenna and collect each focus point.

assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'shape_estimate:MissingConfig', ...
    'Run build_cfg.m before shape_estimate.m.');

cfg.opts.numIterations = 1; % Use one TR cycle per antenna estimate.

scriptDir = fileparts(mfilename('fullpath'));
figsDir = fullfile(scriptDir, 'figs', 'shape_estimate');
if ~isfolder(figsDir)
    mkdir(figsDir);
end

epsrTolerance = 10*eps(max(1, max(abs(cfg.grid.epsr(:)))));
conductivityTolerance = 10*eps(max(1, max(abs(cfg.grid.cond_e(:)))));
targetMask = ...
    abs(cfg.grid.epsr - cfg.grid.background.epsr) > epsrTolerance | ...
    abs(cfg.grid.cond_e - cfg.grid.background.cond_e) > conductivityTolerance;

points = zeros(cfg.antennas.numAntennas, 2);
totalRuntime = 0;

for antennaIdx = 1:cfg.antennas.numAntennas
    runTimer = tic;
    cfg.antennas.txAntennas = antennaIdx;
    cfg.source.samples(:) = 0;
    cfg.source.samples(antennaIdx, :) = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);
    [itrResult, focusPoint] = tr_mwi(cfg);
    points(antennaIdx, :) = focusPoint;
    elapsedTime = toc(runTimer);
    totalRuntime = totalRuntime + elapsedTime;

    saveMaxEz(itrResult, focusPoint, antennaIdx, figsDir, targetMask, cfg);
    fprintf('\nAntenna %d/%d done transmitting.\n', ...
        antennaIdx, cfg.antennas.numAntennas);
end

averageRuntime = totalRuntime/cfg.antennas.numAntennas;
fprintf('Total runtime: %.2f seconds.\n', totalRuntime);
fprintf('Average runtime per antenna: %.2f seconds.\n', averageRuntime);


%% Plot all antenna focus estimates
saveFocusOverview(points, cfg, targetMask, figsDir);

%% Plotting functions

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

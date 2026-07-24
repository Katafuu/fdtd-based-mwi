function [pointsFigure, outlineFigure] = plotShapeEstimateOverview( ...
        points, cfg, targetMask, sigmaThreshold)
%plotShapeEstimateOverview Plot point and outline shape-estimate diagnostics.

arguments
    points (:, 2) double
    cfg (1, 1) struct
    targetMask (:, :) logical
    sigmaThreshold (1, 1) double {mustBeFinite, mustBePositive} = 1.5
end

validateattributes(targetMask, {'logical'}, {'size', [cfg.Nx cfg.Ny]}, ...
    mfilename, 'targetMask');
validateattributes(cfg.antennas.doiMask, {'logical'}, ...
    {'size', [cfg.Nx cfg.Ny]}, mfilename, 'cfg.antennas.doiMask');

[doiPoints, doiPointIndices, outsideDoiPoints] = ...
    selectShapePointsInsideDoi(points, cfg.antennas.doiMask);
[inlierPoints, outlierPoints, inlierDoiIndices] = ...
    filterShapeOutlinePoints(doiPoints, sigmaThreshold);
inlierOriginalIndices = doiPointIndices(inlierDoiIndices);
[orderedPoints, ~] = ...
    orderShapeOutlinePoints(inlierPoints, inlierOriginalIndices);

[pointsFigure, pointsAxes] = createBaseFigure( ...
    cfg, targetMask, 'shapeEstimatePoints');
hold(pointsAxes, 'on');

if ~isempty(inlierPoints)
    inlierHandle = plot(pointsAxes, ...
        inlierPoints(:, 1), inlierPoints(:, 2), 'wo', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.5);
    set(inlierHandle, ...
        'DisplayName', 'Maximum points', 'Tag', 'inlierPoints');
end

if ~isempty(outlierPoints)
    outlierHandle = plot(pointsAxes, ...
        outlierPoints(:, 1), outlierPoints(:, 2), 'x', ...
        'Color', [0.65 0.65 0.65], 'MarkerSize', 6, 'LineWidth', 1);
    set(outlierHandle, ...
        'DisplayName', 'Excluded statistical outliers', 'Tag', 'outlierPoints');
end

hold(pointsAxes, 'off');
title(pointsAxes, sprintf( ...
    'Focus points (%d retained, %d statistical outliers, %d outside DOI)', ...
    size(inlierPoints, 1), size(outlierPoints, 1), size(outsideDoiPoints, 1)));
legend(pointsAxes, 'show', 'Location', 'bestoutside');

[outlineFigure, outlineAxes] = createBaseFigure( ...
    cfg, targetMask, 'shapeEstimateOutline');
hold(outlineAxes, 'on');

if size(orderedPoints, 1) >= 3
    closedPoints = [orderedPoints; orderedPoints(1, :)];
    estimateHandle = plot(outlineAxes, ...
        closedPoints(:, 1), closedPoints(:, 2), '-', ...
        'Color', [0 0.9 0.9], 'LineWidth', 1.75);
    set(estimateHandle, ...
        'DisplayName', 'Estimated outline', 'Tag', 'estimatedOutline');
else
    warning('plotShapeEstimateOverview:InsufficientInliers', ...
        'At least three inlier points are required to draw an estimated outline.');
end

hold(outlineAxes, 'off');
title(outlineAxes, sprintf( ...
    'Estimated outline (%d retained points)', size(inlierPoints, 1)));
legend(outlineAxes, 'show', 'Location', 'bestoutside');
end

function [figureHandle, axesHandle] = createBaseFigure(cfg, targetMask, figureTag)
%createBaseFigure Create one shape-estimate figure with shared context.

figureHandle = figure('Visible', 'off', 'Tag', figureTag);
axesHandle = axes(figureHandle);
set(axesHandle, 'Color', [0.12 0.12 0.12]);
hold(axesHandle, 'on');

[~, doiHandle] = contour(axesHandle, ...
    double(cfg.antennas.doiMask).', [0.5 0.5], ...
    'Color', [0 0.45 0.74], 'LineStyle', '--', 'LineWidth', 1.5);
set(doiHandle, 'DisplayName', 'DOI', 'Tag', 'doiOutline');

[~, targetHandle] = contour(axesHandle, double(targetMask).', [0.5 0.5], ...
    'Color', [1 0.25 0.25], 'LineWidth', 1.75);
set(targetHandle, 'DisplayName', 'True target outline', 'Tag', 'targetOutline');

antennaHandle = plot(axesHandle, ...
    cfg.antennas.pos(:, 1), cfg.antennas.pos(:, 2), ...
    'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
set(antennaHandle, 'DisplayName', 'Antennas', 'Tag', 'antennaPositions');

hold(axesHandle, 'off');
axis(axesHandle, 'equal');
xlim(axesHandle, [1 cfg.Nx]);
ylim(axesHandle, [1 cfg.Ny]);
set(axesHandle, 'YDir', 'normal');
grid(axesHandle, 'on');
xlabel(axesHandle, 'x grid index');
ylabel(axesHandle, 'y grid index');
end

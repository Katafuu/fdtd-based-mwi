function handles = plotDomainOverlays(axesHandle, cfg, txAntenna, lineColor)
%plotDomainOverlays Draw array, transmitter, DOI, and target geometry.

antennaPositions = cfg.antennas.pos;
doiMask = cfg.antennas.doiMask ~= 0;
wasHeld = ishold(axesHandle);
hold(axesHandle, 'on');

handles = struct();
handles.antennas = plot(axesHandle, antennaPositions(:, 1), ...
    antennaPositions(:, 2), 'o', 'Color', lineColor, ...
    'MarkerSize', 6, 'LineWidth', 1.2, ...
    'DisplayName', 'Antennas', 'Tag', 'antennaPositions');

txPosition = antennaPositions(txAntenna, :);
handles.transmitter = plot(axesHandle, txPosition(:, 1), ...
    txPosition(:, 2), 'rx', 'MarkerSize', 11, 'LineWidth', 2.2, ...
    'DisplayName', 'Active TX', 'Tag', 'txAntenna');
[~, handles.doi] = contour(axesHandle, double(doiMask'), [0.5 0.5], ...
    'c--', 'LineWidth', 1.25, 'HandleVisibility', 'off', ...
    'Tag', 'doiOutline');

handles.targets = gobjects(0, 1);
if isfield(cfg, 'targets') && ~isempty(cfg.targets)
    handles.targets = gobjects(numel(cfg.targets), 1);
    for targetIndex = 1:numel(cfg.targets)
        handles.targets(targetIndex) = drawTargetOutline( ...
            axesHandle, cfg.targets(targetIndex), lineColor);
    end
end

if ~wasHeld
    hold(axesHandle, 'off');
end
end

function targetHandle = drawTargetOutline(axesHandle, target, lineColor)
switch string(target.name)
    case "circle"
        theta = linspace(0, 2*pi, 256);
        targetHandle = plot(axesHandle, ...
            target.properties.center(1) + target.properties.radius*cos(theta), ...
            target.properties.center(2) + target.properties.radius*sin(theta), ...
            '-', 'Color', lineColor, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off', 'Tag', 'targetOutline');
    case "rectangle"
        bounds = target.properties.bounds;
        x = [bounds(1), bounds(2), bounds(2), bounds(1), bounds(1)];
        y = [bounds(3), bounds(3), bounds(4), bounds(4), bounds(3)];
        targetHandle = plot(axesHandle, x, y, '-', 'Color', lineColor, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off', ...
            'Tag', 'targetOutline');
    case "triangle"
        vertices = [target.properties.vertices; ...
            target.properties.vertices(1, :)];
        targetHandle = plot(axesHandle, vertices(:, 1), vertices(:, 2), ...
            '-', 'Color', lineColor, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off', 'Tag', 'targetOutline');
    otherwise
        error('plotDomainOverlays:UnsupportedTargetShape', ...
            'Unsupported target shape "%s".', target.name);
end
end

function coverage = plotCircularFootprintCoverage( ...
        mask4D, xRange, yRange, pairAccepted, cfg, focusPoint)
%plotCircularFootprintCoverage Plot accepted TR-shifted mask coverage.
%
% coverage(x,y) counts all accepted directed TX-RX masks. Reciprocal
% TX-RX masks are combined only when drawing their footprint outlines.

validateattributes(mask4D, {'logical'}, {'nonempty'}, ...
    mfilename, 'mask4D');
validateattributes(xRange, {'numeric'}, {'real', 'finite', 'vector'}, ...
    mfilename, 'xRange');
validateattributes(yRange, {'numeric'}, {'real', 'finite', 'vector'}, ...
    mfilename, 'yRange');
validateattributes(pairAccepted, {'logical'}, {'2d', 'nonempty'}, ...
    mfilename, 'pairAccepted');
validateattributes(focusPoint, {'numeric'}, ...
    {'real', 'finite', 'vector', 'numel', 2}, mfilename, 'focusPoint');

numTransmitters = size(mask4D, 3);
numReceivers = size(mask4D, 4);
if size(mask4D, 1) ~= numel(xRange) || size(mask4D, 2) ~= numel(yRange)
    error('plotCircularFootprintCoverage:RangeSizeMismatch', ...
        'xRange and yRange must match the first two dimensions of mask4D.');
end
if ~isequal(size(pairAccepted), [numTransmitters numReceivers])
    error('plotCircularFootprintCoverage:PairSizeMismatch', ...
        'pairAccepted must match the TX-RX dimensions of mask4D.');
end
if numTransmitters ~= numReceivers
    error('plotCircularFootprintCoverage:NonSquarePairGrid', ...
        'Reciprocal outline plotting requires equal TX and RX counts.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'antennas') || ...
        ~isfield(cfg.antennas, 'pos') || ...
        ~isequal(size(cfg.antennas.pos), [numTransmitters 2])
    error('plotCircularFootprintCoverage:InvalidAntennaPositions', ...
        'cfg.antennas.pos must contain one [x y] row per antenna.');
end

acceptedMasks = mask4D & reshape(pairAccepted, ...
    1, 1, numTransmitters, numReceivers);
interiorCoverage = sum(acceptedMasks, [3 4]);
coverage = zeros(cfg.Nx, cfg.Ny);
coverage(xRange, yRange) = interiorCoverage;

axesHandle = gca;
imagesc(axesHandle, 1:cfg.Nx, 1:cfg.Ny, coverage.', ...
    'Tag', 'maskCoverage');
axis(axesHandle, 'equal');
xlim(axesHandle, [0.5 cfg.Nx + 0.5]);
ylim(axesHandle, [0.5 cfg.Ny + 0.5]);
set(axesHandle, 'YDir', 'normal');
xlabel(axesHandle, 'x index');
ylabel(axesHandle, 'y index');
title(axesHandle, 'Accepted TR-shifted footprint coverage');
colormap(axesHandle, turbo);
colorbarHandle = colorbar(axesHandle);
colorbarHandle.Label.String = 'Accepted directed-mask overlap count';

maximumCoverage = max(coverage(:));
if maximumCoverage > 0
    clim(axesHandle, [0 maximumCoverage]);
end

wasHeld = ishold(axesHandle);
hold(axesHandle, 'on');
for tx = 1:numTransmitters-1
    for rx = tx+1:numReceivers
        reciprocalMask = false(size(mask4D, 1), size(mask4D, 2));
        if pairAccepted(tx, rx)
            reciprocalMask = reciprocalMask | mask4D(:, :, tx, rx);
        end
        if pairAccepted(rx, tx)
            reciprocalMask = reciprocalMask | mask4D(:, :, rx, tx);
        end
        if any(reciprocalMask(:))
            contour(axesHandle, xRange, yRange, double(reciprocalMask.'), ...
                [0.5 0.5], 'LineColor', [0.15 0.15 0.15], ...
                'LineWidth', 0.5, 'HandleVisibility', 'off', ...
                'Tag', 'maskOutline');
        end
    end
end

overlayHandles = plotDomainOverlays(axesHandle, cfg, ...
    cfg.antennas.txAntennas(1), 'w');
focusHandle = plot(axesHandle, focusPoint(1), focusPoint(2), 'wp', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', ...
    'MarkerSize', 12, 'LineWidth', 1.5, ...
    'DisplayName', 'TR focus', 'Tag', 'trFocus');
legend(axesHandle, [overlayHandles.antennas ...
    overlayHandles.transmitter focusHandle], 'Location', 'best');

if ~wasHeld
    hold(axesHandle, 'off');
end
end
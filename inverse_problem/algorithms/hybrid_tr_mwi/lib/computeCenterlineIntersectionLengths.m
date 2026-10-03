function intersectionLengths = computeCenterlineIntersectionLengths( ...
        targetMask, antennaPositions, cellSize, pairEligible, sampleDensity)
%computeCenterlineIntersectionLengths Measure target chord lengths per pair.

validateattributes(targetMask, {'logical'}, {'2d', 'nonempty'}, ...
    mfilename, 'targetMask');
validateattributes(antennaPositions, {'numeric'}, ...
    {'real', 'finite', '2d', 'ncols', 2}, mfilename, 'antennaPositions');
cellSize = double(cellSize(:).');
if numel(cellSize) ~= 2 || any(~isfinite(cellSize)) || any(cellSize <= 0)
    error('computeCenterlineIntersectionLengths:InvalidCellSize', ...
        'cellSize must contain two positive finite spacings [dx dy].');
end
numAntennas = size(antennaPositions, 1);
validateattributes(pairEligible, {'logical'}, ...
    {'2d', 'size', [numAntennas numAntennas]}, ...
    mfilename, 'pairEligible');
if nargin < 5 || isempty(sampleDensity)
    sampleDensity = 20;
end
validateattributes(sampleDensity, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'sampleDensity');

gridSize = size(targetMask);
if any(antennaPositions(:, 1) < 1 | ...
        antennaPositions(:, 1) > gridSize(1) | ...
        antennaPositions(:, 2) < 1 | ...
        antennaPositions(:, 2) > gridSize(2))
    error('computeCenterlineIntersectionLengths:AntennaOutOfBounds', ...
        'Every antenna position must lie inside targetMask.');
end

intersectionLengths = nan(numAntennas, numAntennas);
for tx = 1:numAntennas
    for rx = 1:numAntennas
        if tx == rx || ~pairEligible(tx, rx)
            continue
        end

        startPoint = antennaPositions(tx, :);
        endPoint = antennaPositions(rx, :);
        directionCells = endPoint - startPoint;
        segmentLengthCells = hypot(directionCells(1), directionCells(2));
        segmentLengthPhysical = hypot( ...
            directionCells(1)*cellSize(1), ...
            directionCells(2)*cellSize(2));
        numSamples = max(2, ceil(segmentLengthCells*sampleDensity)+1);
        sampleParameter = linspace(0, 1, numSamples);
        sampleX = startPoint(1) + sampleParameter*directionCells(1);
        sampleY = startPoint(2) + sampleParameter*directionCells(2);
        insideTarget = interp2(1:gridSize(1), 1:gridSize(2), ...
            double(targetMask.'), sampleX, sampleY, 'nearest', 0) > 0.5;
        targetFraction = trapz(sampleParameter, double(insideTarget));
        intersectionLengths(tx, rx) = ...
            targetFraction * segmentLengthPhysical;
    end
end
end

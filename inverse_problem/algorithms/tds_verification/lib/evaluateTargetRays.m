function [rayEstimates, targetPixelCounts, backgroundPixelCounts, ...
        pairAccepted, coverage] = evaluateTargetRays( ...
        averageEpsr, backgroundEpsr, antennaPositions, targetMask, ...
        rayWidthPixels, pairEligible)
%evaluateTargetRays Estimate homogeneous target epsr for eligible thick rays.

validateattributes(antennaPositions, {'numeric'}, ...
    {'real', 'finite', '2d', 'ncols', 2}, mfilename, 'antennaPositions');
numAntennas = size(antennaPositions, 1);
validateattributes(averageEpsr, {'numeric'}, ...
    {'real', '2d', 'size', [numAntennas numAntennas]}, ...
    mfilename, 'averageEpsr');
validateattributes(backgroundEpsr, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'backgroundEpsr');
validateattributes(targetMask, {'logical'}, {'2d', 'nonempty'}, ...
    mfilename, 'targetMask');
validateattributes(rayWidthPixels, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'rayWidthPixels');
validateattributes(pairEligible, {'logical'}, ...
    {'2d', 'size', [numAntennas numAntennas]}, ...
    mfilename, 'pairEligible');

gridSize = size(targetMask);
rayEstimates = nan(numAntennas, numAntennas);
targetPixelCounts = zeros(numAntennas, numAntennas);
backgroundPixelCounts = zeros(numAntennas, numAntennas);
pairAccepted = false(numAntennas, numAntennas);
coverage = zeros(gridSize);

for tx = 1:numAntennas
    for rx = 1:numAntennas
        if tx == rx || ~pairEligible(tx, rx)
            continue
        end

        rayMask = buildThickRayMask(gridSize, ...
            antennaPositions(tx, :), antennaPositions(rx, :), ...
            rayWidthPixels);
        numTargetPixels = nnz(rayMask & targetMask);
        if numTargetPixels == 0
            continue
        end

        numBackgroundPixels = nnz(rayMask) - numTargetPixels;
        targetPixelCounts(tx, rx) = numTargetPixels;
        backgroundPixelCounts(tx, rx) = numBackgroundPixels;

        if ~isfinite(averageEpsr(tx, rx))
            continue
        end

        estimate = estimateHomogeneousEpsr(averageEpsr(tx, rx), ...
            backgroundEpsr, numTargetPixels, numBackgroundPixels);
        if ~isfinite(estimate)
            continue
        end

        rayEstimates(tx, rx) = estimate;
        pairAccepted(tx, rx) = true;
        coverage = coverage + double(rayMask);
    end
end
end

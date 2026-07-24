function [inlierPoints, outlierPoints, inlierIndices, outlierIndices, statistics] = ...
        filterShapeOutlinePoints(points, sigmaThreshold)
%filterShapeOutlinePoints Exclude coordinate outliers from shape estimates.

arguments
    points (:, 2) double
    sigmaThreshold (1, 1) double {mustBeFinite, mustBePositive} = 1.5
end

numPoints = size(points, 1);
finiteMask = all(isfinite(points), 2);
inlierMask = false(numPoints, 1);

if any(finiteMask)
    finitePoints = points(finiteMask, :);
    coordinateMean = mean(finitePoints, 1);
    coordinateStd = std(finitePoints, 0, 1);

    finiteInlierMask = true(size(finitePoints, 1), 1);
    varyingCoordinates = coordinateStd > 0;
    if any(varyingCoordinates)
        coordinateDeviation = abs( ...
            finitePoints(:, varyingCoordinates) - coordinateMean(varyingCoordinates));
        coordinateLimit = sigmaThreshold*coordinateStd(varyingCoordinates);
        finiteInlierMask = all(coordinateDeviation <= coordinateLimit, 2);
    end

    inlierMask(finiteMask) = finiteInlierMask;
else
    coordinateMean = [NaN NaN];
    coordinateStd = [NaN NaN];
end

inlierIndices = find(inlierMask);
outlierIndices = find(~inlierMask);
inlierPoints = points(inlierIndices, :);
outlierPoints = points(outlierIndices, :);
statistics = struct( ...
    'Mean', coordinateMean, ...
    'StandardDeviation', coordinateStd, ...
    'SigmaThreshold', sigmaThreshold);
end

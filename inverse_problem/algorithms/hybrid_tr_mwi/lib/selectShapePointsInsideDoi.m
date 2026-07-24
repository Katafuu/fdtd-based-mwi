function [insidePoints, insideIndices, outsidePoints, outsideIndices] = ...
        selectShapePointsInsideDoi(points, doiMask)
%selectShapePointsInsideDoi Keep finite integer points inside the DOI mask.

arguments
    points (:, 2) double
    doiMask (:, :) logical
end

numPoints = size(points, 1);
validIndexMask = all(isfinite(points), 2) & ...
    all(points == round(points), 2) & ...
    points(:, 1) >= 1 & points(:, 1) <= size(doiMask, 1) & ...
    points(:, 2) >= 1 & points(:, 2) <= size(doiMask, 2);

insideMask = false(numPoints, 1);
validIndices = find(validIndexMask);
if ~isempty(validIndices)
    linearIndices = sub2ind(size(doiMask), ...
        points(validIndices, 1), points(validIndices, 2));
    insideMask(validIndices) = doiMask(linearIndices);
end

insideIndices = find(insideMask);
outsideIndices = find(~insideMask);
insidePoints = points(insideIndices, :);
outsidePoints = points(outsideIndices, :);
end

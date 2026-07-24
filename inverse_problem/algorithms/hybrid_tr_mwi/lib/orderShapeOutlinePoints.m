function [orderedPoints, orderedOriginalIndices] = ...
        orderShapeOutlinePoints(points, originalIndices)
%orderShapeOutlinePoints Build a deterministic nearest-unvisited point tour.

arguments
    points (:, 2) double {mustBeFinite}
    originalIndices (:, 1) double {mustBeFinite, mustBePositive, mustBeInteger}
end

numPoints = size(points, 1);
if numel(originalIndices) ~= numPoints
    error('orderShapeOutlinePoints:SizeMismatch', ...
        'originalIndices must contain one index for every point.');
end
if numel(unique(originalIndices)) ~= numPoints
    error('orderShapeOutlinePoints:DuplicateIndices', ...
        'originalIndices must contain unique values.');
end
if numPoints == 0
    orderedPoints = zeros(0, 2);
    orderedOriginalIndices = zeros(0, 1);
    return
end

[~, startOrder] = sortrows([points originalIndices], [1 2 3]);
currentIndex = startOrder(1);
visited = false(numPoints, 1);
tour = zeros(numPoints, 1);
tour(1) = currentIndex;
visited(currentIndex) = true;

for tourPosition = 2:numPoints
    candidateIndices = find(~visited);
    offsets = points(candidateIndices, :) - points(currentIndex, :);
    distanceSquared = sum(offsets.^2, 2);
    minimumDistance = min(distanceSquared);
    tieTolerance = 16*eps(max(1, minimumDistance));
    tiedCandidates = candidateIndices( ...
        abs(distanceSquared - minimumDistance) <= tieTolerance);
    [~, tieOrder] = min(originalIndices(tiedCandidates));
    currentIndex = tiedCandidates(tieOrder);
    tour(tourPosition) = currentIndex;
    visited(currentIndex) = true;
end

orderedPoints = points(tour, :);
orderedOriginalIndices = originalIndices(tour);
end

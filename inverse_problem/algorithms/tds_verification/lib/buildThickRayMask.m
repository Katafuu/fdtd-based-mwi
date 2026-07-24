function rayMask = buildThickRayMask(gridSize, startPoint, endPoint, rayWidthPixels)
%buildThickRayMask Build a finite thick strip between two grid points.

gridSize = double(gridSize(:).');
if numel(gridSize) ~= 2 || any(~isfinite(gridSize)) || ...
        any(gridSize < 1) || any(gridSize ~= round(gridSize))
    error('buildThickRayMask:InvalidGridSize', ...
        'gridSize must contain two positive integers.');
end

startPoint = validatePoint(startPoint, gridSize, ...
    'buildThickRayMask:InvalidStartPoint');
endPoint = validatePoint(endPoint, gridSize, ...
    'buildThickRayMask:InvalidEndPoint');
validateattributes(rayWidthPixels, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'rayWidthPixels');

direction = endPoint - startPoint;
segmentLengthSquared = sum(direction.^2);
if segmentLengthSquared == 0
    error('buildThickRayMask:DegenerateSegment', ...
        'startPoint and endPoint must be different.');
end

[xGrid, yGrid] = ndgrid(1:gridSize(1), 1:gridSize(2));
projection = ((xGrid - startPoint(1)) .* direction(1) + ...
    (yGrid - startPoint(2)) .* direction(2)) ./ segmentLengthSquared;
projection = max(0, min(1, projection));

closestX = startPoint(1) + projection .* direction(1);
closestY = startPoint(2) + projection .* direction(2);
distanceToSegment = hypot(xGrid - closestX, yGrid - closestY);

rayMask = distanceToSegment <= rayWidthPixels/2;
end

function point = validatePoint(point, gridSize, errorId)
point = double(point(:).');
if numel(point) ~= 2 || any(~isfinite(point)) || ...
        point(1) < 1 || point(1) > gridSize(1) || ...
        point(2) < 1 || point(2) > gridSize(2)
    error(errorId, ...
        'The point must be a finite [x y] coordinate inside the grid.');
end
end

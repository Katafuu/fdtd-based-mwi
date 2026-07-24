function region = shape_polygon(grid, vertices, units)
%shape_polygon Create a polygonal region using inpolygon and edge distances.

if nargin < 3 || isempty(units)
    units = "index";
end
if ~isnumeric(vertices) || size(vertices, 2) ~= 2 || size(vertices, 1) < 3 || ...
        any(~isfinite(vertices(:)))
    error('fdtdgeom:shape_polygon:InvalidVertices', 'vertices must be an N-by-2 finite numeric array.');
end

[x, y] = coordinatesForUnits(grid, units);
vertices = closeVertices(double(vertices));

inside = inpolygon(x, y, vertices(:, 1), vertices(:, 2));
distance = distanceToEdges(x, y, vertices);
signedDistance = distance;
signedDistance(inside) = -distance(inside);

region = fdtdgeom.utility_regionFromSignedDistance(signedDistance, "polygon", units);
region.mask = inside;
region.depthInside = max(0, -region.signedDistance);
region.depthOutside = max(0, region.signedDistance);
end

function [x, y] = coordinatesForUnits(grid, units)
units = char(string(units));
if strcmpi(units, 'index')
    x = grid.xIndex;
    y = grid.yIndex;
elseif strcmpi(units, 'physical')
    x = grid.xPhysical;
    y = grid.yPhysical;
else
    error('fdtdgeom:shape_polygon:InvalidUnits', 'units must be "index" or "physical".');
end
end

function vertices = closeVertices(vertices)
if any(vertices(1, :) ~= vertices(end, :))
    vertices(end + 1, :) = vertices(1, :);
end
end

function distance = distanceToEdges(x, y, vertices)
distance = inf(size(x));
for edgeIndex = 1:size(vertices, 1) - 1
    x1 = vertices(edgeIndex, 1);
    y1 = vertices(edgeIndex, 2);
    x2 = vertices(edgeIndex + 1, 1);
    y2 = vertices(edgeIndex + 1, 2);

    vx = x2 - x1;
    vy = y2 - y1;
    segmentLengthSquared = vx * vx + vy * vy;
    t = ((x - x1) * vx + (y - y1) * vy) / segmentLengthSquared;
    t = min(max(t, 0), 1);
    projectionX = x1 + t * vx;
    projectionY = y1 + t * vy;
    distance = min(distance, hypot(x - projectionX, y - projectionY));
end
end

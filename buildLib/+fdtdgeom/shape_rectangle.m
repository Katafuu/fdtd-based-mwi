function region = shape_rectangle(grid, bounds, units)
%shape_rectangle Create an axis-aligned rectangular region.

if nargin < 3 || isempty(units)
    units = "index";
end

bounds = rowVector(bounds, 4, 'fdtdgeom:shape_rectangle:InvalidBounds');
if bounds(1) > bounds(2) || bounds(3) > bounds(4)
    error('fdtdgeom:shape_rectangle:InvalidBounds', 'bounds must be [xmin xmax ymin ymax].');
end

[x, y] = coordinatesForUnits(grid, units);
centerX = 0.5 * (bounds(1) + bounds(2));
centerY = 0.5 * (bounds(3) + bounds(4));
halfWidth = 0.5 * (bounds(2) - bounds(1));
halfHeight = 0.5 * (bounds(4) - bounds(3));

qx = abs(x - centerX) - halfWidth;
qy = abs(y - centerY) - halfHeight;
outsideDistance = hypot(max(qx, 0), max(qy, 0));
insideDistance = min(max(qx, qy), 0);
signedDistance = outsideDistance + insideDistance;

region = fdtdgeom.utility_regionFromSignedDistance(signedDistance, "rectangle", units);
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
    error('fdtdgeom:shape_rectangle:InvalidUnits', 'units must be "index" or "physical".');
end
end

function value = rowVector(value, expectedLength, errorId)
value = double(value(:).');
if numel(value) ~= expectedLength || any(~isfinite(value))
    error(errorId, 'Expected a finite vector with %d elements.', expectedLength);
end
end

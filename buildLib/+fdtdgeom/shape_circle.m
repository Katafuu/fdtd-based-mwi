function region = shape_circle(grid, center, radius, units)
%shape_circle Create a circular region as a mask and signed distance field.

if nargin < 4 || isempty(units)
    units = "index";
end
if ~isscalar(radius) || ~isnumeric(radius) || radius <= 0 || ~isfinite(radius)
    error('fdtdgeom:shape_circle:InvalidRadius', 'radius must be a positive finite scalar.');
end

[x, y] = coordinatesForUnits(grid, units);
center = rowVector(center, 2, 'fdtdgeom:shape_circle:InvalidCenter');

signedDistance = hypot(x - center(1), y - center(2)) - radius;
region = fdtdgeom.utility_regionFromSignedDistance(signedDistance, "circle", units);
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
    error('fdtdgeom:shape_circle:InvalidUnits', 'units must be "index" or "physical".');
end
end

function value = rowVector(value, expectedLength, errorId)
value = double(value(:).');
if numel(value) ~= expectedLength || any(~isfinite(value))
    error(errorId, 'Expected a finite vector with %d elements.', expectedLength);
end
end

function region = shape_annulus(grid, center, innerRadius, outerRadius, units)
%shape_annulus Create a ring-shaped region.

if nargin < 5 || isempty(units)
    units = "index";
end
if ~isscalar(innerRadius) || ~isscalar(outerRadius) || ...
        ~isnumeric(innerRadius) || ~isnumeric(outerRadius) || ...
        innerRadius < 0 || outerRadius <= innerRadius || ...
        ~isfinite(innerRadius) || ~isfinite(outerRadius)
    error('fdtdgeom:shape_annulus:InvalidRadii', 'Require 0 <= innerRadius < outerRadius.');
end

[x, y] = coordinatesForUnits(grid, units);
center = rowVector(center, 2, 'fdtdgeom:shape_annulus:InvalidCenter');

distanceFromCenter = hypot(x - center(1), y - center(2));
signedDistance = max(innerRadius - distanceFromCenter, distanceFromCenter - outerRadius);
region = fdtdgeom.utility_regionFromSignedDistance(signedDistance, "annulus", units);
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
    error('fdtdgeom:shape_annulus:InvalidUnits', 'units must be "index" or "physical".');
end
end

function value = rowVector(value, expectedLength, errorId)
value = double(value(:).');
if numel(value) ~= expectedLength || any(~isfinite(value))
    error(errorId, 'Expected a finite vector with %d elements.', expectedLength);
end
end

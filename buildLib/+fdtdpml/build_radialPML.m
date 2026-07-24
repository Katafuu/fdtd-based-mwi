function pml = build_radialPML(grid, center, thickness, order, sigmaMax, kappaMax, units)
%build_radialPML Create a circular radial PML projected onto Cartesian component maps.

if nargin < 7 || isempty(units)
    units = "index";
end
if ~isscalar(thickness) || thickness <= 0 || ~isfinite(thickness)
    error('fdtdpml:build_radialPML:InvalidThickness', 'thickness must be a positive finite scalar.');
end
if ~isscalar(sigmaMax) || sigmaMax < 0 || ~isfinite(sigmaMax)
    error('fdtdpml:build_radialPML:InvalidSigmaMax', 'sigmaMax must be a nonnegative finite scalar.');
end
if ~isscalar(kappaMax) || kappaMax < 1 || ~isfinite(kappaMax)
    error('fdtdpml:build_radialPML:InvalidKappaMax', 'kappaMax must be a finite scalar greater than or equal to 1.');
end

center = rowVector(center, 2, 'fdtdpml:build_radialPML:InvalidCenter');
[x, y] = coordinatesForUnits(grid, units);

outerRadius = outerRadiusForGrid(x, y, center);
innerRadius = outerRadius - thickness;

if innerRadius < 0
    error('fdtdpml:build_radialPML:InvalidThickness', ...
        'thickness must be smaller than or equal to the inferred outer radius.');
end

outerBoundary = fdtdgeom.shape_circle(grid, center, outerRadius, units);
innerBoundary = fdtdgeom.shape_circle(grid, center, innerRadius, units);

annulusMask = innerBoundary.depthOutside > 0 & outerBoundary.mask;
depth = zeros(size(x));
depth(annulusMask) = min(innerBoundary.depthOutside(annulusMask), thickness);
profile = fdtdpml.utility_buildNormalizedPowerProfile(depth, thickness, order);

deltaX = x - center(1);
deltaY = y - center(2);
radialDistance = hypot(deltaX, deltaY);
cosTheta = zeros(size(radialDistance));
sinTheta = zeros(size(radialDistance));
hasDirection = radialDistance > 0;
cosTheta(hasDirection) = deltaX(hasDirection) ./ radialDistance(hasDirection);
sinTheta(hasDirection) = deltaY(hasDirection) ./ radialDistance(hasDirection);

cosSquared = cosTheta .^ 2;
sinSquared = sinTheta .^ 2;
condR = sigmaMax * profile;
kappaR = 1 + (kappaMax - 1) * profile;

pml = struct();
pml.center = center;
pml.outerRadius = outerRadius;
pml.innerRadius = innerRadius;
pml.thickness = thickness;
pml.units = string(units);
pml.backgroundMask = outerBoundary.mask;
pml.innerMask = innerBoundary.mask;
pml.depth = depth;
pml.profile = profile;
pml.radialDistance = radialDistance;
pml.cosTheta = cosTheta;
pml.sinTheta = sinTheta;
pml.condr = condR;
pml.kappar = kappaR;
pml.condx = condR .* cosSquared;
pml.condy = condR .* sinSquared;
pml.kx = 1 + (kappaR - 1) .* cosSquared;
pml.ky = 1 + (kappaR - 1) .* sinSquared;
pml.mask = profile > 0;
end

function outerRadius = outerRadiusForGrid(x, y, center)
xMin = min(x(:));
xMax = max(x(:));
yMin = min(y(:));
yMax = max(y(:));

axisDistances = [
    center(1) - xMin
    xMax - center(1)
    center(2) - yMin
    yMax - center(2)
];

if any(axisDistances < 0)
    error('fdtdpml:build_radialPML:CenterOutsideGrid', ...
        'center must be inside the grid coordinate bounds.');
end

outerRadius = min(axisDistances);

if outerRadius <= 0 || ~isfinite(outerRadius)
    error('fdtdpml:build_radialPML:InvalidOuterRadius', ...
        'inferred outer radius must be positive and finite.');
end
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
    error('fdtdpml:build_radialPML:InvalidUnits', 'units must be "index" or "physical".');
end
end

function value = rowVector(value, expectedLength, errorId)
value = double(value(:).');
if numel(value) ~= expectedLength || any(~isfinite(value))
    error(errorId, 'Expected a finite vector with %d elements.', expectedLength);
end
end
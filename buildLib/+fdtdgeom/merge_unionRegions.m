function region = merge_unionRegions(varargin)
%merge_unionRegions Combine regions with a logical OR mask.

validateRegionCount(nargin, 'fdtdgeom:merge_unionRegions:NoRegions');
region = varargin{1};
mask = region.mask;
signedDistance = region.signedDistance;

for regionIndex = 2:nargin
    nextRegion = varargin{regionIndex};
    validateCompatible(region, nextRegion, 'fdtdgeom:merge_unionRegions:IncompatibleRegions');
    mask = mask | nextRegion.mask;
    signedDistance = min(signedDistance, nextRegion.signedDistance);
end

region = buildRegion(mask, signedDistance, "union", region.units);
end

function validateRegionCount(count, errorId)
if count < 1
    error(errorId, 'At least one region is required.');
end
end

function validateCompatible(firstRegion, secondRegion, errorId)
if ~isequal(size(firstRegion.mask), size(secondRegion.mask))
    error(errorId, 'All regions must have the same mask size.');
end
end

function region = buildRegion(mask, signedDistance, name, units)
region = struct();
region.mask = mask;
region.signedDistance = signedDistance;
region.depthInside = max(0, -signedDistance);
region.depthOutside = max(0, signedDistance);
region.name = string(name);
region.units = string(units);
end

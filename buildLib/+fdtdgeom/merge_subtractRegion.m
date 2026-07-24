function region = merge_subtractRegion(baseRegion, cutoutRegion)
%merge_subtractRegion Return baseRegion with cutoutRegion removed.

if ~isequal(size(baseRegion.mask), size(cutoutRegion.mask))
    error('fdtdgeom:merge_subtractRegion:IncompatibleRegions', 'Both regions must have the same mask size.');
end

mask = baseRegion.mask & ~cutoutRegion.mask;
signedDistance = max(baseRegion.signedDistance, -cutoutRegion.signedDistance);

region = struct();
region.mask = mask;
region.signedDistance = signedDistance;
region.depthInside = max(0, -signedDistance);
region.depthOutside = max(0, signedDistance);
region.name = "subtraction";
region.units = string(baseRegion.units);
end

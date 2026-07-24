function region = utility_regionFromSignedDistance(signedDistance, name, units)
%utility_regionFromSignedDistance Build a region struct from a signed distance field.

if nargin < 2 || isempty(name)
    name = "";
end
if nargin < 3 || isempty(units)
    units = "index";
end
if ~isnumeric(signedDistance) || ~ismatrix(signedDistance) || any(~isfinite(signedDistance(:)))
    error('fdtdgeom:utility_regionFromSignedDistance:InvalidSignedDistance', ...
        'signedDistance must be a finite numeric matrix.');
end

region = struct();
region.mask = signedDistance <= 0;
region.signedDistance = signedDistance;
region.depthInside = max(0, -signedDistance);
region.depthOutside = max(0, signedDistance);
region.name = string(name);
region.units = string(units);
end

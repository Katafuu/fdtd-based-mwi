function pml = buildPMLProfile(depthX, depthY, thickness, order, sigmaMax, kappaMax)
%buildPMLProfile Build PML conductivity and kappa maps from component depths.

if ~isequal(size(depthX), size(depthY))
    error('fdtdpml:buildPMLProfile:SizeMismatch', 'depthX and depthY must have the same size.');
end
if ~isscalar(sigmaMax) || sigmaMax < 0 || ~isfinite(sigmaMax)
    error('fdtdpml:buildPMLProfile:InvalidSigmaMax', 'sigmaMax must be a nonnegative finite scalar.');
end
if ~isscalar(kappaMax) || kappaMax < 1 || ~isfinite(kappaMax)
    error('fdtdpml:buildPMLProfile:InvalidKappaMax', 'kappaMax must be a finite scalar greater than or equal to 1.');
end

profileX = fdtdpml.utility_buildNormalizedPowerProfile(depthX, thickness, order);
profileY = fdtdpml.utility_buildNormalizedPowerProfile(depthY, thickness, order);

pml = struct();
pml.depthX = depthX;
pml.depthY = depthY;
pml.profileX = profileX;
pml.profileY = profileY;
pml.condx = sigmaMax * profileX;
pml.condy = sigmaMax * profileY;
pml.kx = 1 + (kappaMax - 1) * profileX;
pml.ky = 1 + (kappaMax - 1) * profileY;
pml.mask = profileX > 0 | profileY > 0;
end

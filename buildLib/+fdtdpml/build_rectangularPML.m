function pml = build_rectangularPML(grid, thickness, order, sigmaMax, kappaMax)
%build_rectangularPML Create the Cartesian rectangular PML used by the C TMz solver.

if ~isscalar(thickness) || thickness <= 0 || fix(thickness) ~= thickness
    error('fdtdpml:build_rectangularPML:InvalidThickness', 'thickness must be a positive integer scalar.');
end

sizeXY = grid.sizeXY;
x0 = grid.xIndex - 1;
y0 = grid.yIndex - 1;

distX = min(x0, sizeXY(1) - 1 - x0);
distY = min(y0, sizeXY(2) - 1 - y0);
depthX = max(0, thickness - distX);
depthY = max(0, thickness - distY);

pml = fdtdpml.buildPMLProfile(depthX, depthY, thickness, order, sigmaMax, kappaMax);
end

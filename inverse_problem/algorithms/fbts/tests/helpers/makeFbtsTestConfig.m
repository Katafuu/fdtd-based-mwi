function cfg = makeFbtsTestConfig()
cfg = struct('c0',3e8,'mu0',4*pi*1e-7,'eps0',8.854e-12, ...
    'Nx',48,'Ny',48,'sizeZ',1,'dx',5e-3,'dy',5e-3);
cfg.dt = 1/(cfg.c0*sqrt(1/cfg.dx^2+1/cfg.dy^2));
cfg.deltaF = 1/(200*cfg.dt);
cfg.Nt = ceil(1/(cfg.dt*cfg.deltaF));
cfg.snapshotStart = 0; cfg.snapshotStride = 1;
cfg.returnEz = true; cfg.returnHx = false; cfg.returnHy = false;
cfg.returnRxSignals = true;
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny],[cfg.dx cfg.dy]);
cfg.pml = fdtdpml.build_rectangularPML(cfg.grid,8,3,2,6);
cfg.pml.type = 'cpml'; cfg.pml.enabled = true;
cfg.pml.ax = 1; cfg.pml.ay = 1; cfg.pml.az = 1;
cfg.pml.thickness = 8;
cfg.antennas = struct('numAntennas',4,'txAntennas',1:4, ...
    'center',[24 24],'radius',14,'pmlPadding',2,'focusPadding',4);
cfg = buildCircularAntennaArrayIdx(cfg);
[x,y] = ndgrid(1:cfg.Nx,1:cfg.Ny);
cfg.grid.epsr(hypot(x-23,y-25)<=4) = 2;
cfg.targets = struct('name','circle','properties',struct('center',[23 25], ...
    'radius',4),'material',struct('epsr',2,'cond_e',0));
cfg.source = struct('location',cfg.antennas.pos, ...
    'samples',zeros(cfg.antennas.numAntennas,cfg.Nt));
end

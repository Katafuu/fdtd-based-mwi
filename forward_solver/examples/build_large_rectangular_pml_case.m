% build_fdtd_case Build the shared TMz FDTD case in memory.
%
% This script creates the scalar simulation parameters, material maps, PML
% maps, and a MEX-ready cfg struct. It intentionally performs no file I/O and
% does not run either solver front end.

exampleDir = fileparts(mfilename('fullpath'));
if isempty(exampleDir)
    exampleDir = pwd;
end
solverRoot = fileparts(exampleDir);
workspaceRoot = fileparts(solverRoot);
addpath(fullfile(workspaceRoot, 'buildLib'));

resultsDir = fullfile(solverRoot, 'results');
inputDir = fullfile(resultsDir, 'fdtd_input');

assert(isfile(fullfile(solverRoot, 'Makefile')) && ...
       isfolder(fullfile(solverRoot, 'fdtd')) && ...
       isfolder(fullfile(solverRoot, 'utility')), ...
       'Keep this script under forward_solver/examples.');

eps0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c0 = 1/sqrt(mu0*eps0);

Nx = 450;
Ny = 450;
dx = 0.005; % sweet spot for freq of 1.5GHz to be represented by 20 points per wavelength, aka physical domain bigger, less fine sampling of wave
dy = dx;
maxTime = 900;
sourceFreq = 1.5e9;
dt = 0.99/(c0*sqrt(1/dx^2 + 1/dy^2));

ax = 1.0;
ay = 1.0;
az = 1.0;

background = struct("epsr", 1.0, "murx", 1.0, "mury", 1.0, ...
    "cond_e", 0.0, "cond_m", 0.0);

grid = fdtdmat.createGrid([Nx Ny], [dx dy], background, [0 0]);

highGrid = fdtdmat.createGrid(2 .* [Nx Ny], [dx/2 dy/2], background, [0 0]);

targetCenter = [(Nx - 1) * dx / 2, (Ny - 1) * dy / 2];
targetRadius = 90 * dx;
target = fdtdgeom.shape_circle(highGrid, targetCenter, targetRadius, "physical");
highGrid = fdtdmat.applyRegion(highGrid, target, struct("epsr", 1.0, "cond_e", 0.0));

ezRows = 1:2:(2*Nx - 1);
ezCols = 1:2:(2*Ny - 1);
hxRows = 1:2:(2*Nx - 1);
hxCols = 2:2:(2*Ny - 2);
hyRows = 2:2:(2*Nx - 2);
hyCols = 1:2:(2*Ny - 1);

grid.epsr = highGrid.epsr(ezRows, ezCols);
grid.cond_e = highGrid.cond_e(ezRows, ezCols);
grid.cond_m = highGrid.cond_m(ezRows, ezCols);
grid.murx = highGrid.murx(hxRows, hxCols);
grid.mury = highGrid.mury(hyRows, hyCols);

pmlThickness = 80;
pmlOrder = 6;
targetReflection = 1e-60;
kappaMax = 6.0;
etaBackground = sqrt((mu0 * grid.background.murx) / (eps0 * grid.background.epsr));
pmlCenter = [(Nx + 1)/2 (Ny + 1)/2];

pmlPhysicalThickness = pmlThickness*dx;
sigmaMax = fdtdpml.utility_sigmaMaxFromReflection(pmlOrder, ...
    targetReflection, etaBackground, pmlPhysicalThickness);
% pml = fdtdpml.build_radialPML(grid, pmlCenter, pmlThickness, ...
    % pmlOrder, sigmaMax, kappaMax, "index");
pml = fdtdpml.build_rectangularPML(grid, pmlThickness, pmlOrder, sigmaMax, kappaMax);

cfg = struct();
cfg.Nx = Nx;
cfg.Ny = Ny;
cfg.sizeZ = 1;
cfg.Nt = maxTime;
cfg.dx = dx;
cfg.dy = dy;
cfg.dt = dt;
cfg.snapshotStart = 0; % which 
cfg.snapshotStride = 1;
cfg.returnEz = true;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = false;

cfg.grid = grid;

cfg.pml = pml;
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = ax;
cfg.pml.ay = ay;
cfg.pml.az = az;
cfg.pml.thickness = pmlThickness;
cfg.pml.m = pmlOrder;
cfg.pml.R = targetReflection;

cfg.antennas = struct();
cfg.antennas.numAntennas = 23;
cfg.antennas.pmlPadding = 5;
cfg.antennas.focusPadding = 20;
cfg.antennas.center = [round(cfg.Nx/2), round(cfg.Ny/2)];
cfg.antennas.radius = floor(min(cfg.grid.sizeXY)/2) - ...
    cfg.pml.thickness - cfg.antennas.pmlPadding;
cfg = buildCircularAntennaArrayIdx(cfg);

cfg.source = struct();
cfg.source.frequency = sourceFreq;
cfg.source.location = cfg.antennas.pos;
cfg.source.func = @(physicalTime) ...
    (1 - 2 .* (pi .* (sourceFreq .* physicalTime - 1)).^2) .* ...
    exp(-(pi .* (sourceFreq .* physicalTime - 1)).^2);
sourceWaveform = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);
cfg.source.samples = repmat(sourceWaveform, cfg.antennas.numAntennas, 1);

disp("Built FDTD case in workspace: cfg, grid, pml, Nx, Ny, dx, dy, dt.");

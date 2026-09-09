%[text] # Build TMz FDTD Case
%[text] This Live Script builds the shared simulation case used by both front ends. It creates scalar parameters, material maps, PML maps, antennas, sampled sources, and the MEX-ready |cfg| struct, but it does not write files or run a solver.
%[text] The workspace outputs are the canonical |cfg| plus normal-orientation |grid| and |pml| construction objects.
exampleDir = fileparts(mfilename('fullpath'));
if isempty(exampleDir)
    exampleDir = pwd;
end
solverRoot = fileparts("C:\Users\Galax\Desktop\PolitoCode");
workspaceRoot = fileparts(solverRoot);
addpath(fullfile(workspaceRoot, 'buildLib')); %[output:6ef51303]
resultsDir = fullfile(solverRoot, 'results');
inputDir = fullfile(resultsDir, 'fdtd_input');
assert(isfile(fullfile(solverRoot, 'Makefile')) && ...
       isfolder(fullfile(solverRoot, 'fdtd')) && ...
       isfolder(fullfile(solverRoot, 'utility')), ...
       'Keep this script under forward_solver/examples.');
%%
%[text] ## Initialize Configuration
%[text] Initialize the canonical configuration before constructing the grid, materials, CPML maps, source, and runtime settings.
cfg = struct();
%%
%[text] ## Basic Grid Parameters
%[text] Define the original free-space constants, grid dimensions, cell spacing, duration, and Courant-limited time step directly in |cfg|.
cfg.eps0 = 8.854187817e-12;
cfg.mu0 = 4*pi*1e-7;
cfg.c0 = 1/sqrt(cfg.mu0*cfg.eps0);
cfg.Nx = 200;
cfg.Ny = 200;
cfg.sizeZ = 1;
cfg.dx = 0.001;
cfg.dy = cfg.dx;
cfg.Nt = 4300;
cfg.dt = 0.99/(cfg.c0*sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
%%
%[text] ## Background Material
%[text] Construct the canonical normal-orientation MATLAB grid retained directly in |cfg.grid|.
background = struct("epsr", 1.0, "murx", 1.0, "mury", 1.0, ...
    "cond_e", 0.0, "cond_m", 0.0);
grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy], background, [0 0]);
%%
%[text] ## High-Resolution Geometry
%[text] Build geometry on the original 2x refined grid, then sample it onto the staggered Yee locations without changing the construction.
highGrid = fdtdmat.createGrid(2 .* [cfg.Nx cfg.Ny], ...
    [cfg.dx/2 cfg.dy/2], background, [0 0]);
targetCenter = [(cfg.Nx - 1) * cfg.dx / 2, (cfg.Ny - 1) * cfg.dy / 2];
targetRadius = 90 * cfg.dx;
target = fdtdgeom.shape_circle(highGrid, targetCenter, targetRadius, "physical");
highGrid = fdtdmat.applyRegion(highGrid, target, struct("epsr", 1.0, "cond_e", 0.0));
%%
%[text] ## Yee-Grid Projection
%[text] Project the refined material maps onto the original Ez, Hx, and Hy storage locations.
ezRows = 1:2:(2*cfg.Nx - 1);
ezCols = 1:2:(2*cfg.Ny - 1);
hxRows = 1:2:(2*cfg.Nx - 1);
hxCols = 2:2:(2*cfg.Ny - 2);
hyRows = 2:2:(2*cfg.Nx - 2);
hyCols = 1:2:(2*cfg.Ny - 1);
grid.epsr = highGrid.epsr(ezRows, ezCols);
grid.cond_e = highGrid.cond_e(ezRows, ezCols);
grid.cond_m = highGrid.cond_m(ezRows, ezCols);
grid.murx = highGrid.murx(hxRows, hxCols);
grid.mury = highGrid.mury(hyRows, hyCols);


% What happened to smoothing?
%%
%[text] ## CPML Profiles
%[text] Preserve the radial CPML construction, including all three a parameters, while grouping its canonical maps and metadata under |cfg.pml|.
cfg.pml = struct();
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = 1.0;
cfg.pml.ay = 1.0;
cfg.pml.az = 1.0;
cfg.pml.thickness = 20;
cfg.pml.m = 3;
cfg.pml.R = 1e-60;
cfg.pml.kappaMax = 4.0;
cfg.pml.etaBackground = sqrt((cfg.mu0 * grid.background.murx) / ...
    (cfg.eps0 * grid.background.epsr));
cfg.pml.pmlCenter = [(cfg.Nx + 1)/2 (cfg.Ny + 1)/2];
cfg.pml.pmlPhysicalThickness = cfg.pml.thickness*cfg.dx;
cfg.pml.sigma_max = fdtdpml.utility_sigmaMaxFromReflection(cfg.pml.m, ...
    cfg.pml.R, cfg.pml.etaBackground, cfg.pml.pmlPhysicalThickness);
pml = fdtdpml.build_radialPML(grid, cfg.pml.pmlCenter, cfg.pml.thickness, ...
    cfg.pml.m, cfg.pml.sigma_max, cfg.pml.kappaMax, "index");
%%
%[text] ## Antenna Array
%[text] Build the same one-based circular antenna description used by the MATLAB solvers.
cfg.grid = grid;
cfg.antennas = struct();
cfg.antennas.numAntennas = 23;
cfg.antennas.pmlPadding = 5;
cfg.antennas.focusPadding = 20;
cfg.antennas.center = [round(cfg.Nx/2), round(cfg.Ny/2)];
cfg.antennas.radius = floor(min(cfg.grid.sizeXY)/2) - ...
    cfg.pml.thickness - cfg.antennas.pmlPadding;
cfg = buildCircularAntennaArrayIdx(cfg);
%%
%[text] ## Sampled Source
%[text] Sample the historical 1.5 GHz Ricker pulse in MATLAB and activate it at every antenna.
sourceFrequency = 1.5e9;
cfg.source = struct();
cfg.source.frequency = sourceFrequency;
cfg.source.location = cfg.antennas.pos;
cfg.source.func = @(physicalTime) ...
    (1 - 2 .* (pi .* (sourceFrequency .* physicalTime - 1)).^2) .* ...
    exp(-(pi .* (sourceFrequency .* physicalTime - 1)).^2);
sourceWaveform = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);
cfg.source.samples = repmat(sourceWaveform, cfg.antennas.numAntennas, 1);
%%
%[text] ## Runtime and Outputs
%[text] Preserve the original final-field MEX snapshot behavior and requested output components.
cfg.snapshotStart = 0;
cfg.snapshotStride = 0;
cfg.returnEz = true;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = false;
%%
%[text] ## Canonical Grid and CPML Maps
%[text] The MEX gateway accepts the normal MATLAB \[x,y\] orientation used by the shared builders.
cfg.grid = grid;
cfg.pml.condx = pml.condx;
cfg.pml.condy = pml.condy;
cfg.pml.kx = pml.kx;
cfg.pml.ky = pml.ky;
%%
%[text] ## Workspace Result
%[text] The case is ready for either the MEX or standalone file runner. The normal |grid| and |pml| locals are retained for the file-backed interface.
disp("Built FDTD case in workspace: cfg, grid, pml."); %[output:5367d3b7]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:6ef51303]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Name is nonexistent or not a directory: \/home\/aly\/Desktop\/fdtd-based-mwi\/forward_solver\/buildLib"}}
%---
%[output:5367d3b7]
%   data: {"dataType":"text","outputData":{"text":"Built FDTD case in workspace: cfg, grid, pml.\n","truncated":false}}
%---

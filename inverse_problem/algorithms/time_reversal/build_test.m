%% Build centered-target 2-D TMz FDTD test case
% This script builds a compact test problem for tr.mlx.
% It leaves cfg in the workspace and does not run the FDTD solver.

algorithmDir = fileparts(mfilename('fullpath'));
if isempty(algorithmDir)
    algorithmDir = pwd;
end
workspaceRoot = fileparts(fileparts(fileparts(algorithmDir)));
addpath(fullfile(workspaceRoot, 'buildLib'));
addpath(fullfile(algorithmDir, 'lib'));

%% 0. Initialize cfg
cfg = struct();

%% 1. Basic grid parameters
cfg.c0 = 3e8;
cfg.mu0 = 4*pi*1e-7;
cfg.eps0 = 8.854e-12;
cfg.Nx = 400;
cfg.Ny = 400;
cfg.dx = 1e-3;
cfg.dy = 1e-3;
cfg.dt = 1 / (cfg.c0 * sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
cfg.Nt = 800;
cfg.sizeZ = 1;
cfg.snapshotStart = 0;
cfg.snapshotStride = 0;
cfg.returnEz = false;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = true;

%% 2. Grid
cfg.grid = struct();
cfg.grid.background = struct( ...
    'epsr', 1.0, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.0, ...
    'cond_m', 0.0);
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy], ...
    cfg.grid.background, [0 0]);

%% 3. PML
cfg.pml = fdtdpml.build_rectangularPML(cfg.grid, 20, 3, ...
    -(3 + 1) * log(1e-60) / ...
    (2 * sqrt(cfg.mu0/cfg.eps0) * 20 * cfg.dx), 1.0);
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = 1.0;
cfg.pml.ay = 1.0;
cfg.pml.az = 1.0;
cfg.pml.thickness = 40;
cfg.pml.m = 3;
cfg.pml.R = 1e-60;
cfg.pml.eta0 = sqrt(cfg.mu0/cfg.eps0);
cfg.pml.sigma_max = -(cfg.pml.m + 1) * log(cfg.pml.R) / ...
    (2 * cfg.pml.eta0 * cfg.pml.thickness * cfg.dx);
cfg.pml.sigma_e = max(cfg.pml.condx, cfg.pml.condy);

%% 4. Targets
cfg.targets = struct();
cfg.targets.name = 'circle';
cfg.targets.properties = struct( ...
    'center', [round(cfg.Nx/2), round(cfg.Ny/2)], ...
    'radius', 20);
cfg.targets.material = struct( ...
    'epsr', 5.0, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.1, ...
    'cond_m', 0.0);
cfg.targets.mask = fdtdgeom.shape_circle(cfg.grid, ...
    cfg.targets.properties.center, cfg.targets.properties.radius, 'index').mask;
cfg.grid = fdtdmat.applyRegion(cfg.grid, cfg.targets, cfg.targets.material);

%% 5. Antennas
cfg.antennas.numAntennas = 12;
cfg.antennas.txAntennas = 12; % Select the initial TR transmitting antennas.
cfg.antennas.pmlPadding = 5;
cfg.antennas.focusPadding = 10;
cfg.antennas.orientation = 'top';
cfg = buildPlanarAntennaIdx(cfg);

%% 6. Source
pulse_width = 20e-12;
t_delay = 4 * pulse_width;
amplitude = 10;
t = (0:cfg.Nt-1) .* cfg.dt;
source_normalization = max(abs(-((t - t_delay) ./ pulse_width.^2) .* ...
    exp(-0.5 .* ((t - t_delay) ./ pulse_width).^2)));
cfg.source.location = cfg.antennas.pos;
cfg.source.func = @(physicalTime) amplitude .* ...
    (-((physicalTime - t_delay) ./ pulse_width.^2) .* ...
    exp(-0.5 .* ((physicalTime - t_delay) ./ pulse_width).^2)) ./ source_normalization;
cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);

cfg.opts = struct();
cfg.opts.storeFieldHistory = false; % Skip optional field-history retention.
cfg.filename = '';

disp('Built centered-target FDTD test case in workspace. Run tr.mlx next.');
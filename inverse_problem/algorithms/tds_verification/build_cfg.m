% build_cfg Construct the reproducible transmission-delay verification case.
% Run this script once before tds_2d.m.

%% Repository / path setup
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end

algorithmDir = scriptDir;
algorithmsDir = fileparts(algorithmDir);
workspaceRoot = fileparts(fileparts(fileparts(algorithmDir)));

buildLibDir = fullfile(workspaceRoot, 'buildLib');
mexDir = fullfile(workspaceRoot, 'forward_solver', 'mex');
mwiLibDir = fullfile(algorithmsDir, 'transmission_mwi', 'lib');
verificationLibDir = fullfile(algorithmDir, 'lib');

addpath(buildLibDir, '-end');
addpath(mexDir, '-end');
addpath(mwiLibDir, '-end');
addpath(verificationLibDir, '-begin');
addpath(algorithmDir, '-begin');

assert(isfile(fullfile(workspaceRoot, 'forward_solver', 'Makefile')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'fdtd')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'utility')), ...
       'Could not find forward_solver at the workspace root.');

%% Reproducible homogeneous target controls
if ~exist('verificationSeed', 'var') || isempty(verificationSeed)
    verificationSeed = 1;
end
randomSeed = verificationSeed;
rng(randomSeed, 'twister');

targetOpts = struct();
targetOpts.numTargetsRange = [1 1];
targetOpts.allowedShapes = "rectangle";
targetOpts.radiusRange = [10 28];
targetOpts.sideRange = [31 31];
targetOpts.epsrRange = [2.5 6.0];
targetOpts.condRange = [0.02 0.15];
targetOpts.maxAttempts = 800;
targetOpts.allowOverlap = false;

%% Simulation parameters
cfg = struct();
cfg.c0 = 3e8;
cfg.mu0 = 4*pi*1e-7;
cfg.eps0 = 8.854e-12;
cfg.Nx = 400;
cfg.Ny = 400;
cfg.sizeZ = 1;
cfg.Nt = 600;
cfg.dx = 1e-3;
cfg.dy = 1e-3;
cfg.dt = 1 / (cfg.c0 * sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
cfg.snapshotStart = 0;
cfg.snapshotStride = 1;
cfg.returnEz = true;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = true;

%% Grid
cfg.grid = struct();
cfg.grid.background = struct( ...
    'epsr', 1.0, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.0, ...
    'cond_m', 0.0);
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy], ...
    cfg.grid.background, [0 0]);
highGrid = fdtdmat.createGrid(2 .* [cfg.Nx cfg.Ny], ...
    [cfg.dx/2 cfg.dy/2], cfg.grid.background, [0 0]);

%% PML construction
cfg.pml = fdtdpml.build_rectangularPML(cfg.grid, floor(cfg.Ny * 0.25), 6, ...
    -(6 + 1) * log(1e-60) / ...
    (2 * sqrt((cfg.mu0 * cfg.grid.background.murx) / ...
    (cfg.eps0 * cfg.grid.background.epsr)) * floor(cfg.Ny * 0.25) * cfg.dx), 6.0);
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = 1.0;
cfg.pml.ay = 1.0;
cfg.pml.az = 1.0;
cfg.pml.thicknessRatio = 0.25;
cfg.pml.thickness = floor(cfg.Ny * cfg.pml.thicknessRatio);
cfg.pml.m = 6;
cfg.pml.R = 1e-60;
cfg.pml.kappaMax = 6.0;
cfg.pml.eta0 = sqrt((cfg.mu0 * cfg.grid.background.murx) / ...
    (cfg.eps0 * cfg.grid.background.epsr));
cfg.pml.physicalThickness = cfg.pml.thickness * cfg.dx;
cfg.pml.sigma_max = -(cfg.pml.m + 1) * log(cfg.pml.R) / ...
    (2 * cfg.pml.eta0 * cfg.pml.physicalThickness);
cfg.pml.sigma_e = max(cfg.pml.condx, cfg.pml.condy);

%% Circular antenna array
cfg.antennas = struct();
cfg.antennas.numAntennas = 16;
cfg.antennas.txAntennas = 12;
cfg.antennas.pmlPadding = 5;
cfg.antennas.focusPadding = 20;
cfg.antennas.center = [round(cfg.Nx/2), round(cfg.Ny/2)];
cfg.antennas.xMin = cfg.pml.thickness + 1 + cfg.antennas.pmlPadding;
cfg.antennas.xMax = cfg.Nx - cfg.pml.thickness - cfg.antennas.pmlPadding;
cfg.antennas.yMin = cfg.pml.thickness + 1 + cfg.antennas.pmlPadding;
cfg.antennas.yMax = cfg.Ny - cfg.pml.thickness - cfg.antennas.pmlPadding;
cfg.antennas.radius = floor(min([ ...
    cfg.antennas.center(1) - cfg.antennas.xMin, ...
    cfg.antennas.xMax - cfg.antennas.center(1), ...
    cfg.antennas.center(2) - cfg.antennas.yMin, ...
    cfg.antennas.yMax - cfg.antennas.center(2)]));
cfg = buildCircularAntennaArrayIdx(cfg);

fprintf('Built circular array: Nant = %d, radius = %.1f cells.\n', ...
    cfg.antennas.numAntennas, cfg.antennas.radius);

%% One random homogeneous target on the 2x material grid
cfg.targets = generateRandomTargetSpecs( ...
    cfg.grid, cfg.antennas.doiMask, targetOpts);

assert(isscalar(cfg.targets), ...
    'build_cfg:ExpectedSingleTarget', ...
    'The verification configuration requires exactly one target.');

target = cfg.targets(1);
switch string(target.name)
    case "circle"
        centerPhysical = (target.properties.center - 1) .* ...
            [cfg.dx cfg.dy];
        radiusPhysical = target.properties.radius * cfg.dx;
        radiusPhysical = radiusPhysical + ...
            10 * eps(max([abs(centerPhysical), radiusPhysical]));
        highGridTarget = fdtdgeom.shape_circle(highGrid, ...
            centerPhysical, radiusPhysical, 'physical');
    case "rectangle"
        bounds = target.properties.bounds;
        physicalBounds = [ ...
            (bounds(1) - 1) * cfg.dx, (bounds(2) - 1) * cfg.dx, ...
            (bounds(3) - 1) * cfg.dy, (bounds(4) - 1) * cfg.dy];
        boundsTolerance = 10 * eps(max(abs(physicalBounds)));
        physicalBounds = physicalBounds + ...
            [-boundsTolerance boundsTolerance ...
            -boundsTolerance boundsTolerance];
        highGridTarget = fdtdgeom.shape_rectangle( ...
            highGrid, physicalBounds, 'physical');
    case "triangle"
        physicalVertices = (target.properties.vertices - 1) .* ...
            [cfg.dx cfg.dy];
        highGridTarget = fdtdgeom.shape_polygon( ...
            highGrid, physicalVertices, 'physical');
    otherwise
        error('build_cfg:UnsupportedTargetShape', ...
            'Unsupported target shape "%s".', target.name);
end
highGrid = fdtdmat.applyRegion(highGrid, highGridTarget, target.material);

cfg.grid.epsr = highGrid.epsr(1:2:(2*cfg.Nx - 1), 1:2:(2*cfg.Ny - 1));
cfg.grid.cond_e = highGrid.cond_e(1:2:(2*cfg.Nx - 1), 1:2:(2*cfg.Ny - 1));
cfg.grid.cond_m = highGrid.cond_m(1:2:(2*cfg.Nx - 1), 1:2:(2*cfg.Ny - 1));
cfg.grid.murx = highGrid.murx(1:2:(2*cfg.Nx - 1), 2:2:(2*cfg.Ny - 2));
cfg.grid.mury = highGrid.mury(2:2:(2*cfg.Nx - 2), 1:2:(2*cfg.Ny - 1));

targetName = string({cfg.targets.name}).';
targetEpsr = arrayfun(@(targetSpec) targetSpec.material.epsr, cfg.targets).';
targetCondE = arrayfun(@(targetSpec) targetSpec.material.cond_e, cfg.targets).';
targetTable = table(targetName, targetEpsr, targetCondE, ...
    'VariableNames', {'Name', 'EpsR', 'CondE'});
disp(targetTable);

%% Sampled source
cfg.source = struct();
cfg.source.pulseWidth = 20e-12;
cfg.source.delay = 4 * cfg.source.pulseWidth;
cfg.source.amplitude = 10;
cfg.source.time = (0:cfg.Nt-1) .* cfg.dt;
cfg.source.normalization = max(abs( ...
    -((cfg.source.time - cfg.source.delay) ./ cfg.source.pulseWidth.^2) .* ...
    exp(-0.5 .* ((cfg.source.time - cfg.source.delay) ./ cfg.source.pulseWidth).^2)));
cfg.source.location = cfg.antennas.pos;
cfg.source.func = @(physicalTime) ...
    cfg.source.amplitude .* ...
    (-((physicalTime - cfg.source.delay) ./ cfg.source.pulseWidth.^2) .* ...
    exp(-0.5 .* ((physicalTime - cfg.source.delay) ./ cfg.source.pulseWidth).^2)) ./ ...
    cfg.source.normalization;
cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);
cfg.source.samples(cfg.antennas.txAntennas, :) = repmat( ...
    cfg.source.func(cfg.source.time), ...
    numel(cfg.antennas.txAntennas), 1);

cfg.filename = '';

fprintf('Verification target uses fixed RNG seed %d.\n', randomSeed);

% build_cfg Construct the lossless FBTS FDTD configuration.
% Run this script from main.m before runFbts.m.

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
trLibDir = fullfile(algorithmsDir, 'time_reversal', 'lib');

addpath(buildLibDir, '-end');
addpath(mexDir, '-end');
addpath(trLibDir, '-end');
addpath(algorithmDir, '-begin');

resultsDir = fullfile(algorithmDir, 'results');
inputDir = fullfile(resultsDir, 'mats');

assert(isfile(fullfile(workspaceRoot, 'forward_solver', 'Makefile')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'fdtd')) && ...
       isfolder(fullfile(workspaceRoot, 'forward_solver', 'utility')), ...
       'Could not find forward_solver at the workspace root.');

eps0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c0 = 1/sqrt(mu0*eps0);

%% Simulation parameters
cfg = struct();
cfg.c0 = 3e8;
cfg.mu0 = 4*pi*1e-7;
cfg.eps0 = 8.854e-12;
cfg.Nx = 280;
cfg.Ny = 280;
cfg.sizeZ = 1;
cfg.dx = 1e-3;
cfg.dy = 1e-3;
cfg.dt = 1 / (cfg.c0 * sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
cfg.deltaF = 1e9 / sqrt(2);
cfg.Nt = max(ceil(1 / (cfg.dt * cfg.deltaF)), ...
    ceil(6e-9 / cfg.dt) + 1);
cfg.snapshotStart = 0;
cfg.snapshotStride = 1;
cfg.returnEz = true;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = true;

%% Grid
cfg.grid = struct();
cfg.grid.background = struct( ...
    'epsr', 45.0, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.0, ...
    'cond_m', 0.0);
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy], ...
    cfg.grid.background, [0 0]);
highGrid = fdtdmat.createGrid(2 .* [cfg.Nx cfg.Ny], ...
    [cfg.dx/2 cfg.dy/2], cfg.grid.background, [0 0]);


%% PML construction
pmlThickness = 40;
cfg.pml = fdtdpml.build_rectangularPML(cfg.grid, pmlThickness, 6, ...
    -(6 + 1) * log(1e-60) / ...
    (2 * sqrt((cfg.mu0 * cfg.grid.background.murx) / ...
    (cfg.eps0 * cfg.grid.background.epsr)) * pmlThickness * cfg.dx), 6.0);
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = 1.0;
cfg.pml.ay = 1.0;
cfg.pml.az = 1.0;
cfg.pml.thicknessRatio = pmlThickness / cfg.Ny;
cfg.pml.thickness = pmlThickness;
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
% One shared antenna list. Antenna tx acts as transmitter; all other
% antennas are receivers for that tx.
cfg.antennas = struct();
% Change this count to generate a different circular array and scenario set.
cfg.antennas.numAntennas = 12;
cfg.antennas.txAntennas = 1:cfg.antennas.numAntennas;
cfg.antennas.pmlPadding = 5;
cfg.antennas.focusPadding = 19;
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

%% Available centered targets on the 2x material grid
targetMaterial = struct( ...
    'epsr', 2.0, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.0, ...
    'cond_m', 0.0);
targetCenter = cfg.antennas.center;
circleRadiusCells = 0.01 / cfg.dx;
triangleSidePhysical = 0.017320508075688773;
triangleHalfSideCells = triangleSidePhysical / (2 * cfg.dx);
triangleHeightCells = sqrt(3) * triangleSidePhysical / (2 * cfg.dy);
squareHalfWidthCells = 0.01 / cfg.dx;
squareHalfHeightCells = 0.01 / cfg.dy;

cfg.availableTargets = repmat(struct( ...
    'mask', [], ...
    'name', "", ...
    'material', targetMaterial, ...
    'properties', struct()), 1, 3);
cfg.availableTargets(1).name = "circle";
cfg.availableTargets(1).properties = struct( ...
    'center', targetCenter, 'radius', circleRadiusCells);
cfg.availableTargets(2).name = "triangle";
cfg.availableTargets(2).properties = struct('vertices', [ ...
    targetCenter(1), targetCenter(2) - 2 * triangleHeightCells / 3; ...
    targetCenter(1) - triangleHalfSideCells, ...
        targetCenter(2) + triangleHeightCells / 3; ...
    targetCenter(1) + triangleHalfSideCells, ...
        targetCenter(2) + triangleHeightCells / 3]);
cfg.availableTargets(3).name = "square";
cfg.availableTargets(3).properties = struct('bounds', [ ...
    targetCenter(1) - squareHalfWidthCells, ...
    targetCenter(1) + squareHalfWidthCells, ...
    targetCenter(2) - squareHalfHeightCells, ...
    targetCenter(2) + squareHalfHeightCells]);
cfg.targets = cfg.availableTargets(1);

for targetIdx = 1:numel(cfg.targets)
    target = cfg.targets(targetIdx);
    switch string(target.name)
        case "circle"
            centerPhysical = (target.properties.center - 1) .* ...
                [cfg.dx cfg.dy];
            radiusPhysical = target.properties.radius * cfg.dx;
            radiusPhysical = radiusPhysical + ...
                10 * eps(max([abs(centerPhysical), radiusPhysical]));
            highGridTarget = fdtdgeom.shape_circle(highGrid, ...
                centerPhysical, radiusPhysical, 'physical');
        case "square"
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
    cfg.targets(targetIdx).mask = highGridTarget.mask( ...
        1:2:(2*cfg.Nx - 1), 1:2:(2*cfg.Ny - 1));
    highGrid = fdtdmat.applyRegion( ...
        highGrid, highGridTarget, target.material);
end

cfg.grid.epsr = highGrid.epsr(1:2:(2*cfg.Nx - 1), 1:2:(2*cfg.Ny - 1));
cfg.grid.cond_e = highGrid.cond_e(1:2:(2*cfg.Nx - 1), 1:2:(2*cfg.Ny - 1));
cfg.grid.cond_m = highGrid.cond_m(1:2:(2*cfg.Nx - 1), 1:2:(2*cfg.Ny - 1));

targetName = string({cfg.targets.name}).';
targetEpsr = arrayfun(@(targetSpec) targetSpec.material.epsr, cfg.targets).';
targetCondE = arrayfun(@(targetSpec) targetSpec.material.cond_e, cfg.targets).';
targetTable = table(targetName, targetEpsr, targetCondE, ...
    'VariableNames', {'Name', 'EpsR', 'CondE'});
disp(targetTable);

% The paper pulse is shared by synthetic measurement generation and runFbts.
cfg.source = struct();
cfg.source.location = cfg.antennas.pos;
cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);
tau = 0.125e-9;
cfg.source.func = @(t) ...
    (4 .* t.^3 ./ tau.^4 - t.^4 ./ tau.^5) .* exp(-t ./ tau);

% Algorithm controls are declared directly in runFbts.m.

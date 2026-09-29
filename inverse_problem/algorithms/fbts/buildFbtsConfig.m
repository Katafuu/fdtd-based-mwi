function cfg = buildFbtsConfig(options)
%buildFbtsConfig Reproducible lossless scenes; see README for builder options.
if nargin < 1, options = struct(); end
defaults = struct('seed', [], 'targetOptions', struct(), 'Nx',400,'Ny',400, ...
    'dx',1e-3,'dy',1e-3,'deltaF',1e9/sqrt(2),'numAntennas',8, ...
    'pmlThicknessRatio',0.25,'pmlPadding',5,'focusPadding',20, ...
    'backgroundPermittivity',1,'doiDiameter',[],'targetSpecs',[]);
assert(isstruct(options) && isscalar(options),'fbts:InvalidBuilderOptions','Expected a scalar options struct.');
unknown = setdiff(fieldnames(options),fieldnames(defaults));
assert(isempty(unknown),'fbts:InvalidBuilderOptions','Unknown builder option.');
for key = fieldnames(options).', defaults.(key{1})=options.(key{1}); end
options = defaults;
for key = {'Nx','Ny','numAntennas'}
    validateattributes(options.(key{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
for key = {'dx','dy','deltaF'}
    validateattributes(options.(key{1}),{'numeric'},{'scalar','positive','finite'});
end
validateattributes(options.pmlThicknessRatio,{'numeric'},{'scalar','>',0,'<',0.5});
validateattributes(options.pmlPadding,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(options.focusPadding,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(options.backgroundPermittivity,{'numeric'},{'scalar','finite','>=',1});
if ~isempty(options.doiDiameter)
    validateattributes(options.doiDiameter,{'numeric'},{'scalar','finite','positive'});
end
assert(isstruct(options.targetOptions) && isscalar(options.targetOptions), ...
    'fbts:InvalidBuilderOptions','targetOptions must be a scalar struct.');
if isempty(options.seed), options.seed = randi(2^32-1); end
validateattributes(options.seed,{'numeric'},{'scalar','integer','nonnegative','<=',2^32-1});
previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(options.seed,'twister');
generationState = rng;

% build_cfg Construct the lossless FBTS FDTD configuration.
% Run this script immediately before fbts_demo.m or fbts_batch.m.

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

%% Random target controls
randomSeed = options.seed;
targetOpts = struct();
targetOpts.numTargetsRange = [1 1];
targetOpts.allowedShapes = "rectangle";
targetOpts.radiusRange = [10 28];
targetOpts.sideRange = [31 31];
targetOpts.epsrRange = [2.5 6.0];
targetOpts.condRange = [0 0];
targetOpts.maxAttempts = 800;
targetOpts.allowOverlap = false;
for key = fieldnames(options.targetOptions).'
    assert(isfield(targetOpts,key{1}),'fbts:InvalidBuilderOptions','Unknown target option %s.',key{1});
    targetOpts.(key{1}) = options.targetOptions.(key{1});
end

%% Simulation parameters
cfg = struct();
cfg.c0 = 3e8;
cfg.mu0 = 4*pi*1e-7;
cfg.eps0 = 8.854e-12;
cfg.Nx = options.Nx;
cfg.Ny = options.Ny;
cfg.sizeZ = 1;
cfg.dx = options.dx;
cfg.dy = options.dy;
cfg.dt = 1 / (cfg.c0 * sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
cfg.deltaF = options.deltaF;
cfg.Nt = ceil(1 / (cfg.dt * cfg.deltaF));
cfg.snapshotStart = 0;
cfg.snapshotStride = 1;
cfg.returnEz = true;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = true;

%% Grid
cfg.grid = struct();
cfg.grid.background = struct( ...
    'epsr', options.backgroundPermittivity, ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', 0.0, ...
    'cond_m', 0.0);
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy], ...
    cfg.grid.background, [0 0]);
highGrid = fdtdmat.createGrid(2 .* [cfg.Nx cfg.Ny], ...
    [cfg.dx/2 cfg.dy/2], cfg.grid.background, [0 0]);


%% PML construction
cfg.pml = fdtdpml.build_rectangularPML(cfg.grid, floor(min(cfg.Nx,cfg.Ny) * options.pmlThicknessRatio), 6, ...
    -(6 + 1) * log(1e-60) / ...
    (2 * sqrt((cfg.mu0 * cfg.grid.background.murx) / ...
    (cfg.eps0 * cfg.grid.background.epsr)) * floor(min(cfg.Nx,cfg.Ny) * options.pmlThicknessRatio) * cfg.dx), 6.0);
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = 1.0;
cfg.pml.ay = 1.0;
cfg.pml.az = 1.0;
cfg.pml.thicknessRatio = options.pmlThicknessRatio;
cfg.pml.thickness = floor(min(cfg.Nx,cfg.Ny) * cfg.pml.thicknessRatio);
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
% One shared antenna list. Antenna tx acts as transmitter; all
% antennas, including tx itself, are receivers for that tx.
cfg.antennas = struct();
cfg.antennas.numAntennas = options.numAntennas;
cfg.antennas.txAntennas = 1:cfg.antennas.numAntennas;
cfg.antennas.pmlPadding = options.pmlPadding;
cfg.antennas.focusPadding = options.focusPadding;
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
if ~isempty(options.doiDiameter)
    [ix,iy] = ndgrid(1:cfg.Nx,1:cfg.Ny);
    radius = options.doiDiameter/2;
    assert(radius < cfg.antennas.radius*min(cfg.dx,cfg.dy), ...
        'fbts:InvalidDOI','The DOI must lie strictly inside the antenna ring.');
    cfg.antennas.doiMask = hypot((ix-cfg.antennas.center(1))*cfg.dx, ...
        (iy-cfg.antennas.center(2))*cfg.dy) <= radius;
    cfg.antennas.doiDiameter_m = options.doiDiameter;
    cfg.antennas.focusPadding = cfg.antennas.radius-radius/cfg.dx;
end

fprintf('Built circular array: Nant = %d, radius = %.1f cells.\n', ...
    cfg.antennas.numAntennas, cfg.antennas.radius);

%% Random targets on the 2x material grid
if isempty(options.targetSpecs)
    cfg.targets = generateRandomTargetSpecs(cfg.grid, cfg.antennas.doiMask, targetOpts);
else
    assert(isstruct(options.targetSpecs),'fbts:InvalidTargets','targetSpecs must be a struct array.');
    cfg.targets = options.targetSpecs;
end

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
        case {"rectangle","square"}
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
        case {"triangle","hexagon","polygon"}
            physicalVertices = (target.properties.vertices - 1) .* ...
                [cfg.dx cfg.dy];
            highGridTarget = fdtdgeom.shape_polygon( ...
                highGrid, physicalVertices, 'physical');
        otherwise
            error('build_cfg:UnsupportedTargetShape', ...
                'Unsupported target shape "%s".', target.name);
    end
    cfg.targets(targetIdx).mask = highGridTarget.mask(1:2:(2*cfg.Nx-1),1:2:(2*cfg.Ny-1));
    assert(~any(cfg.targets(targetIdx).mask & ~cfg.antennas.doiMask,'all'), ...
        'fbts:TargetOutsideDOI','The rasterized target must fit inside the DOI.');
    cfg.targets(targetIdx).properties.orientation_deg = 0;
    switch string(target.name)
        case "circle"
            physical = struct('center_m',centerPhysical,'radius_m',radiusPhysical);
        case {"rectangle","square"}
            physical = struct('bounds_m',physicalBounds, ...
                'center_m',[mean(physicalBounds(1:2)),mean(physicalBounds(3:4))]);
        case {"triangle","hexagon","polygon"}
            physical = struct('vertices_m',physicalVertices,'center_m',mean(physicalVertices,1));
    end
    cfg.targets(targetIdx).physical = physical;
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

% runFbts defines the paper pulse and activates one transmitter per solve.
cfg.source = struct();
cfg.source.location = cfg.antennas.pos;
cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);

% Algorithm controls are resolved by resolveFbtsOptions and may be
% overridden through the optional fbtsOptions workspace structure.

cfg.generation = struct('seed',randomSeed,'rng_algorithm','twister', ...
    'rng_before_targets',generationState,'rng_after_targets',rng, ...
    'builder_options',options,'target_options',targetOpts, ...
    'material_grid_factor',2,'sampling','odd-index point sampling', ...
    'target_masks_source','same high-grid rasterization and sampling as material maps');
end

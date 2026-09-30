% main Run the FBTS target, position, and antenna-layout batch.
mainTimer = tic;
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'build_cfg.m'));
run(fullfile(scriptDir, 'gen_combinatorial_scenarios.m'));
baseCfg = cfg;

ensureFdtdMex(workspaceRoot);

batchOutputDir = fullfile(scriptDir, 'batch_output');
if ~isfolder(batchOutputDir)
    [created, message] = mkdir(batchOutputDir);
    if ~created
        error('fbts:CreateBatchOutputFailed', '%s', message);
    end
end

numTargets = numel(baseCfg.availableTargets);
numPositions = size(targetPositions, 1);
numAntennas = baseCfg.antennas.numAntennas;
numLayouts = sum(cellfun(@(rows) size(rows, 1), antennaOffConfigs));
numCases = numTargets * numPositions * numLayouts;
pairIndices = zeros(numCases, 1);
disabledByCase = cell(numCases, 1);
resultFiles = cell(numCases, 1);
reservedPaths = containers.Map('KeyType', 'char', 'ValueType', 'logical');

caseIndex = 0;
for targetIndex = 1:numTargets
    for positionIndex = 1:numPositions
        pairIndex = (targetIndex - 1) * numPositions + positionIndex;
        for offCount = 0:numAntennas-1
            layouts = antennaOffConfigs{offCount + 1};
            for layoutIndex = 1:size(layouts, 1)
                caseIndex = caseIndex + 1;
                disabled = layouts(layoutIndex, :);
                pairIndices(caseIndex) = pairIndex;
                disabledByCase{caseIndex} = disabled;
                if isempty(disabled)
                    disabledLabel = 'all-on';
                else
                    disabledLabel = strjoin(string(disabled), '-');
                end
                baseName = sprintf('fbts_%s_%d_%s', ...
                    char(baseCfg.availableTargets(targetIndex).name), ...
                    positionIndex, char(disabledLabel));
                resultFiles{caseIndex} = nextResultPath( ...
                    batchOutputDir, baseName, reservedPaths);
            end
        end
    end
end

% Size the process pool from available physical memory on the host.
if ispc
    [~, systemMemory] = memory;
    availableBytes = systemMemory.PhysicalMemory.Available;
elseif isunix && isfile('/proc/meminfo')
    memInfo = fileread('/proc/meminfo');
    availableKiB = regexp(memInfo, ...
        'MemAvailable:\s+([0-9]+)\s+kB', 'tokens', 'once');
    if isempty(availableKiB)
        error('fbts:LinuxMemoryUnavailable', ...
            'Could not read MemAvailable from /proc/meminfo.');
    end
    availableBytes = str2double(availableKiB{1}) * 1024;
else
    error('fbts:UnsupportedMemoryPlatform', ...
        'Worker sizing requires Windows or Linux memory information.');
end
bytesPerWorker = 12 * 2^30;
ramWorkerLimit = floor(0.75 * availableBytes / bytesPerWorker);
if ramWorkerLimit < 1
    error('fbts:InsufficientRAM', ...
        'Available physical RAM is below the reserved 12 GiB worker budget.');
end

% One complete transmitter/receiver measurement set per target and position.
pairCfgs = cell(numTargets * numPositions, 1);
pairMeasurements = cell(numTargets * numPositions, 1);
measurementTimer = tic;
for targetIndex = 1:numTargets
    for positionIndex = 1:numPositions
        pairTimer = tic;
        pairIndex = (targetIndex - 1) * numPositions + positionIndex;
        pairCfgs{pairIndex} = configureTarget( ...
            baseCfg, targetIndex, targetPositions(positionIndex, :));
        pairMeasurements{pairIndex} = generateMeasurements(pairCfgs{pairIndex});
        fprintf('Measured target %d/%d, position %d/%d in %.3f s.\n', ...
            targetIndex, numTargets, positionIndex, numPositions, ...
            toc(pairTimer));
    end
end
measurementTime = toc(measurementTimer);

processCluster = parcluster('Processes');
workerCount = min([ramWorkerLimit, processCluster.NumWorkers, numCases]);
existingPool = gcp('nocreate');
if isempty(existingPool)
    parpool('Processes', workerCount);
elseif isa(existingPool, 'parallel.ProcessPool')
    workerCount = min(workerCount, existingPool.NumWorkers);
else
    error('fbts:ProcessPoolRequired', ...
        'The existing parallel pool must use process workers.');
end
fprintf('Using %d process workers for %d FBTS cases.\n', ...
    workerCount, numCases);
fprintf('Shared synthetic measurement time: %.3f seconds.\n', ...
    measurementTime);

progressQueue = parallel.pool.DataQueue;
afterEach(progressQueue, @(info) fprintf( ...
    'Completed case %d/%d (%s): full run %.3f s, MAT write %.3f s, %d bytes.\n', ...
    info.index, numCases, info.file, info.caseTime, ...
    info.writeTime, info.fileBytes));
parfor (caseIndex = 1:numCases, workerCount)
    caseTimer = tic;
    pairIndex = pairIndices(caseIndex);
    disabled = disabledByCase{caseIndex};
    active = setdiff(1:numAntennas, disabled, 'stable');
    caseCfg = pairCfgs{pairIndex};
    fullMeasurements = pairMeasurements{pairIndex};
    EzMeasured = fullMeasurements(active, active, :);

    caseCfg.antennas.pos = caseCfg.antennas.pos(active, :);
    caseCfg.antennas.numAntennas = numel(active);
    caseCfg.antennas.txAntennas = 1:numel(active);
    caseCfg.source.location = caseCfg.antennas.pos;
    caseCfg.source.samples = zeros(numel(active), caseCfg.Nt);

    [results, caseCfg] = runFbts(caseCfg, EzMeasured);
    resultFile = resultFiles{caseIndex};
    results.output_directory = string(batchOutputDir);
    results.output_files = struct('mat_file', string(resultFile));

    writeTimer = tic;
    save_results(resultFile, caseCfg, results);
    writeTime = toc(writeTimer);
    fileInfo = dir(resultFile);
    caseTime = toc(caseTimer);
    send(progressQueue, struct( ...
        'index', caseIndex, 'file', resultFile, ...
        'caseTime', caseTime, 'writeTime', writeTime, ...
        'fileBytes', fileInfo.bytes));
end

fprintf('Whole-main execution time: %.3f seconds.\n', toc(mainTimer));

function ensureFdtdMex(workspaceRoot)
if exist('fdtd_mex', 'file') == 3
    return
end
fprintf('FDTD MEX not found; building it for this platform.\n');
run(fullfile(workspaceRoot, 'forward_solver', 'mex', ...
    'build_fdtd_mex.m'));
rehash;
assert(exist('fdtd_mex', 'file') == 3, ...
    'fbts:MissingFdtdMex', ...
    'The FDTD MEX build completed, but fdtd_mex is not on the MATLAB path.');
end

function caseCfg = configureTarget(baseCfg, targetIndex, targetPosition)
caseCfg = baseCfg;
target = baseCfg.availableTargets(targetIndex);
baseCenterPhysical = (baseCfg.antennas.center - 1) .* ...
    [baseCfg.dx baseCfg.dy];
offsetCells = (targetPosition - baseCenterPhysical) ./ ...
    [baseCfg.dx baseCfg.dy];

highGrid = fdtdmat.createGrid(2 .* [baseCfg.Nx baseCfg.Ny], ...
    [baseCfg.dx/2 baseCfg.dy/2], baseCfg.grid.background, [0 0]);
switch string(target.name)
    case "circle"
        target.properties.center = target.properties.center + offsetCells;
        centerPhysical = (target.properties.center - 1) .* ...
            [baseCfg.dx baseCfg.dy];
        radiusPhysical = target.properties.radius * baseCfg.dx;
        radiusPhysical = radiusPhysical + ...
            10 * eps(max([abs(centerPhysical), radiusPhysical]));
        highGridTarget = fdtdgeom.shape_circle( ...
            highGrid, centerPhysical, radiusPhysical, 'physical');
    case "triangle"
        target.properties.vertices = target.properties.vertices + offsetCells;
        physicalVertices = (target.properties.vertices - 1) .* ...
            [baseCfg.dx baseCfg.dy];
        highGridTarget = fdtdgeom.shape_polygon( ...
            highGrid, physicalVertices, 'physical');
    case "square"
        target.properties.bounds = target.properties.bounds + ...
            [offsetCells(1) offsetCells(1) offsetCells(2) offsetCells(2)];
        bounds = target.properties.bounds;
        physicalBounds = [ ...
            (bounds(1) - 1) * baseCfg.dx, ...
            (bounds(2) - 1) * baseCfg.dx, ...
            (bounds(3) - 1) * baseCfg.dy, ...
            (bounds(4) - 1) * baseCfg.dy];
        boundsTolerance = 10 * eps(max(abs(physicalBounds)));
        physicalBounds = physicalBounds + ...
            [-boundsTolerance boundsTolerance ...
            -boundsTolerance boundsTolerance];
        highGridTarget = fdtdgeom.shape_rectangle( ...
            highGrid, physicalBounds, 'physical');
    otherwise
        error('fbts:UnsupportedTargetShape', ...
            'Unsupported target shape "%s".', target.name);
end

target.mask = highGridTarget.mask( ...
    1:2:(2*baseCfg.Nx - 1), 1:2:(2*baseCfg.Ny - 1));
caseCfg.targets = target;
highGrid = fdtdmat.applyRegion( ...
    highGrid, highGridTarget, target.material);
caseCfg.grid.epsr = highGrid.epsr( ...
    1:2:(2*baseCfg.Nx - 1), 1:2:(2*baseCfg.Ny - 1));
caseCfg.grid.cond_e = highGrid.cond_e( ...
    1:2:(2*baseCfg.Nx - 1), 1:2:(2*baseCfg.Ny - 1));
caseCfg.grid.cond_m = highGrid.cond_m( ...
    1:2:(2*baseCfg.Nx - 1), 1:2:(2*baseCfg.Ny - 1));
end

function EzMeasured = generateMeasurements(cfg)
numAntennas = cfg.antennas.numAntennas;
EzMeasured = zeros(numAntennas, numAntennas, cfg.Nt);
measurementCfg = cfg;
measurementCfg.returnEz = false;
measurementCfg.returnHx = false;
measurementCfg.returnHy = false;
measurementCfg.returnRxSignals = true;
sourcePulse = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);
for transmitter = 1:numAntennas
    measurementCfg.source.samples(:) = 0;
    measurementCfg.source.samples(transmitter, :) = sourcePulse;
    measuredResult = fdtd_mex(measurementCfg);
    EzMeasured(transmitter, :, :) = reshape( ...
        measuredResult.rx_signals, 1, numAntennas, cfg.Nt);
end
end

function resultFile = nextResultPath(outputDir, baseName, reservedPaths)
suffix = 0;
while true
    if suffix == 0
        fileName = [baseName '.mat'];
    else
        fileName = sprintf('%s(%d).mat', baseName, suffix);
    end
    resultFile = fullfile(outputDir, fileName);
    if ~isfile(resultFile) && ~isKey(reservedPaths, resultFile)
        reservedPaths(resultFile) = true;
        return
    end
    suffix = suffix + 1;
end
end

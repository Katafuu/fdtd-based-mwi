function [summary, root, manifest] = runFbtsBatch(scenes, fbtsOptions, batchOptions)
%runFbtsBatch Reflection/rotation classes, shared scans, independent inversions.
%   scenes is one cfg, a struct array, or a cell array of cfg structures.
%   Use batchOptions.dryRun for the complete catalog and resource estimates.
if nargin < 2, fbtsOptions = []; end
if nargin < 3, batchOptions = []; end
preflightTimer = tic;
options = resolveFbtsOptions(fbtsOptions);
batchOptions = resolveFbtsBatchOptions(batchOptions);
if batchOptions.maxWorkers > 1 && ~batchOptions.dryRun
    assert(fbtsParallelAvailable(),'fbts:ParallelUnavailable', ...
        'maxWorkers > 1 requires Parallel Computing Toolbox.');
end
if isstruct(scenes), scenes = num2cell(scenes); end
assert(iscell(scenes) && ~isempty(scenes),'fbts:InvalidScenes','Provide at least one scene configuration.');
scenes = scenes(:).';
for j = 1:numel(scenes)
    scenes{j}.solverThreads = options.solverThreads;
    scenes{j}.solverBackend = options.solverBackend;
    scenes{j}.gpuDeviceIndex = options.gpuDeviceIndex;
end
assert(~strcmp(options.solverBackend,'cuda') || batchOptions.maxWorkers==1, ...
    'fbts:SingleGpuWorker','Use one worker with the single-GPU CUDA backend.');
catalogs = cell(size(scenes)); hashes = cell(size(scenes)); ids = cell(size(scenes));
rows = repmat(fbtsBatchSummaryRow(),0,1);
estimates = repmat(struct('scene_id','','cases',0,'measurement_calls',0, ...
    'fine_iteration_calls',0,'coarse_iteration_calls',0, ...
    'single_Ez_history_bytes',0,'trace_and_image_history_bytes',0),numel(scenes),1);
for j = 1:numel(scenes)
    cfg = scenes{j}; validateFbtsScene(cfg,options);
    ids{j} = sprintf('scene_%06d',j);
    hashes{j} = fbtsHash(fbtsPortableConfig(cfg));
    catalogs{j} = buildFbtsAntennaCatalog(cfg,batchOptions.symmetry, ...
        batchOptions.disabledCounts,batchOptions.seed,ids{j},batchOptions.maxCandidates);
    m = arrayfun(@(entry) numel(entry.active_indices),catalogs{j});
    coarse = prepareCoarseSensitivityCfg(cfg,options.sensitivityDownsampleFactor);
    estimates(j) = struct('scene_id',ids{j},'cases',numel(m), ...
        'measurement_calls',cfg.antennas.numAntennas, ...
        'fine_iteration_calls',2*options.numIterations*sum(m), ...
        'coarse_iteration_calls',2*options.numIterations*sum(m), ...
        'single_Ez_history_bytes',8*cfg.Nx*cfg.Ny*cfg.Nt, ...
        'trace_and_image_history_bytes',8*((cfg.Nt+coarse.Nt)*options.numIterations*sum(m.^2) + ...
        nnz(cfg.antennas.doiMask)*(options.numIterations+1)*numel(m)));
    for k = 1:numel(catalogs{j})
        entry = catalogs{j}(k); row = fbtsBatchSummaryRow();
        row.run_index = numel(rows)+1; row.scene_id = string(ids{j}); row.class_id = string(entry.class_id);
        row.num_active_antennas = numel(entry.active_indices);
        row.active_antenna_indices = join(string(entry.active_indices),';');
        row.disabled_antenna_indices = join(string(entry.disabled_indices),';');
        if isempty(entry.disabled_indices), row.disabled_antenna_indices = ""; end
        row.orbit_size = entry.orbit_size; row.rotation = entry.rotation; row.reflected = entry.reflected;
        row.first_transmitter = entry.active_indices(1); row.num_iterations = options.numIterations;
        row.output_directory = string(fullfile('scenes',ids{j},'cases',entry.class_id));
        row.result_file = string(fullfile(row.output_directory,'result.mat'));
        rows(end+1,1) = row; %#ok<AGROW>
    end
end
numericalOptions = rmfield(options,{'outputDirectory','runLabel'});
selectionOptions = struct('symmetry',batchOptions.symmetry,'seed',batchOptions.seed, ...
    'disabledCounts',batchOptions.disabledCounts);
provenance = fbtsCodeProvenance();
identity = fbtsHash({hashes,numericalOptions,selectionOptions,provenance.code_hash});
manifest = struct('schema_version',1,'identity',identity, ...
    'options',options,'batch_options',batchOptions,'catalogs',{catalogs}, ...
    'scene_hashes',{hashes},'scene_ids',{ids},'scene_references',{cell(size(scenes))}, ...
    'resource_estimates',estimates,'source_code_hash',provenance.code_hash, ...
    'symmetry_caveat','Layout classes only; fixed-scene reconstruction quality need not be invariant.', ...
    'status_rows',rows,'invocation_wall_seconds',[]);
manifest.scene_definitions=cell(size(scenes));
for j=1:numel(scenes)
    definition=struct();
    for name={'study','generation','targets'}
        if isfield(scenes{j},name{1}), definition.(name{1})=scenes{j}.(name{1}); end
    end
    if isfield(definition,'targets') && isfield(definition.targets,'mask')
        definition.targets=rmfield(definition.targets,'mask');
    end
    manifest.scene_definitions{j}=definition;
end
fprintf('FBTS batch: %d scenes, %d configurations, symmetry=%s.\n',numel(scenes),numel(rows),batchOptions.symmetry);
disp(struct2table(estimates));
root = "";
if batchOptions.dryRun, summary = struct2table(rows); return; end
timer = tic;
batchTiming = struct('preflight_seconds',toc(preflightTimer),'worker_startup_seconds',0, ...
    'scene_storage_seconds',zeros(numel(scenes),1),'acquisition_load_or_run_seconds',zeros(numel(scenes),1), ...
    'index_write_seconds',0);
if batchOptions.resume
    assert(strlength(options.outputDirectory)>0 && isfile(fullfile(options.outputDirectory,'manifest.mat')), ...
        'fbts:MissingManifest','Resume requires outputDirectory containing a saved manifest.');
    root = options.outputDirectory;
else
    root = string(createFbtsOutputDirectory(options.outputDirectory, ...
        fullfile(fileparts(mfilename('fullpath')),'figs'),'batch'));
end
root = string(char(java.io.File(char(root)).getCanonicalPath()));
batchLock = acquireFbtsBatchLock(char(root)); %#ok<NASGU>
manifestFile = fullfile(root,'manifest.mat');
if batchOptions.resume
    saved = load(manifestFile,'manifest');
    assert(saved.manifest.schema_version==1 && strcmp(saved.manifest.identity,identity), ...
        'fbts:BatchMismatch','Saved batch differs in scenes, algorithm options, selection, or source/solver code.');
    manifest = saved.manifest;
    % Saved representatives are authoritative; resuming must never redraw them.
    catalogs = manifest.catalogs;
    rows = manifest.status_rows;
    for j = 1:numel(manifest.scene_references)
        refs = manifest.scene_references{j};
        if isempty(refs), continue; end
        for item = {'system','truth','acquisition'}
            assert(strcmp(fbtsHash(fullfile(root,refs.([item{1} '_file'])),'file'), ...
                refs.([item{1} '_hash'])), 'fbts:DatasetMismatch', ...
                'A shared %s file is missing or modified.',item{1});
        end
    end
    % Recover committed statuses even if the coordinator was interrupted before CSV refresh.
    for k = 1:numel(rows)
        path = fullfile(root,rows(k).output_directory,'status.mat');
        if isfile(path)
            stored = load(path,'status');
            assert(stored.status.scene_id==rows(k).scene_id && stored.status.class_id==rows(k).class_id, ...
                'fbts:BatchMismatch','Case status IDs differ.');
            rows(k) = stored.status;
        end
        if rows(k).status=="success"
            assert(isfile(fullfile(root,rows(k).result_file)) && ...
                strcmp(fbtsHash(fullfile(root,rows(k).result_file),'file'),rows(k).result_hash), ...
                'fbts:DatasetMismatch','A completed result is missing or modified.');
            if ~isfile(path)
                saveFbtsAtomic(path,struct('status',rows(k)));
            end
        end
    end
else
    provenance = fbtsCodeProvenance(fullfile(root,'provenance','source'));
    provenance.execution = batchOptions;
    provenance.hardware = fbtsHardwareInfo();
    provenance.precision = 'double';
    provenance.solver_backend = options.solverBackend;
    provenance.solver_threads = options.solverThreads;
    provenance.logical_processors = java.lang.Runtime.getRuntime().availableProcessors();
    provenance.memory_note = 'Measure peak MATLAB/MEX memory before increasing maxWorkers.';
    if isfile('/proc/meminfo')
        info = fileread('/proc/meminfo');
        total = regexp(info,'MemTotal:\s+(\d+) kB','tokens','once');
        if ~isempty(total), provenance.physical_memory_bytes = 1024*str2double(total{1}); end
    end
    saveFbtsAtomic(fullfile(root,'provenance','metadata.mat'),struct('provenance',provenance));
end
sceneRows = cell(numel(scenes),7);
for j = 1:numel(scenes)
    c = scenes{j}; shape = ''; location = NaN; center = [NaN NaN]; sizeM = NaN;
    if isfield(c,'study')
        shape = c.study.shape; location = c.study.location_number; center = c.study.center_m; sizeM = c.study.nominal_size_m;
    end
    sceneRows(j,:) = {ids{j},shape,location,center(1),center(2),sizeM,fullfile('scenes',ids{j},'truth.mat')};
end
writetable(cell2table(sceneRows,'VariableNames',{'scene_id','shape','location_number', ...
    'center_x_m','center_y_m','nominal_size_m','truth_file'}),fullfile(root,'scenes.csv'));
persist();
attempted = 0;
for j = 1:numel(scenes)
    indices = find([rows.scene_id] == string(ids{j}));
    pending = indices([rows(indices).status] ~= "success");
    if isempty(pending) || attempted >= batchOptions.maxCases, continue; end
    sceneTimer = tic;
    references = saveFbtsScene(char(root),ids{j},scenes{j},options);
    batchTiming.scene_storage_seconds(j) = toc(sceneTimer);
    acquisitionTimer = tic;
    acquisitionPath = fullfile(root,references.acquisition_file);
    try
        if isfile(acquisitionPath)
            acquisition = load(acquisitionPath);
            measurements = struct('rx',acquisition.rx_full,'metadata',acquisition.measurement, ...
                'timing',acquisition.measurementTiming);
            selectFbtsMeasurements(measurements,scenes{j});
            previous = manifest.scene_references{j};
            if ~isempty(previous) && isfield(previous,'acquisition_hash')
                assert(strcmp(fbtsHash(acquisitionPath,'file'),previous.acquisition_hash), ...
                    'fbts:DatasetMismatch','Saved acquisition file has been modified.');
            end
        else
            measurements = acquireFbtsMeasurements(scenes{j});
            saveFbtsAtomic(acquisitionPath,struct('rx_full',measurements.rx, ...
                'measurement',measurements.metadata,'measurementTiming',measurements.timing));
        end
        references.acquisition_hash = fbtsHash(acquisitionPath,'file');
    catch exception
        failure = struct('stage','acquisition','identifier',exception.identifier, ...
            'message',exception.message,'stack',exception.stack);
        saveFbtsAtomic(fullfile(root,'scenes',ids{j},'acquisition_failure.mat'),struct('failure',failure));
        for k = pending
            rows(k).status = "failed"; rows(k).failure_stage = "acquisition";
            rows(k).error_id = string(exception.identifier); rows(k).error_message = string(exception.message);
            saveFbtsAtomic(fullfile(root,rows(k).output_directory,'status.mat'),struct('status',rows(k)));
        end
        persist(); continue
    end
    batchTiming.acquisition_load_or_run_seconds(j) = toc(acquisitionTimer);
    manifest.scene_references{j} = references; persist();
    pending = pending(1:min(numel(pending),batchOptions.maxCases-attempted));
    if batchOptions.maxWorkers > 1
        assert(fbtsParallelAvailable(),'fbts:ParallelUnavailable', ...
            'maxWorkers > 1 requires Parallel Computing Toolbox.');
        poolTimer = tic;
        pool = gcp('nocreate');
        if isempty(pool), pool = parpool('Processes',batchOptions.maxWorkers); end
        assert(~isa(pool,'parallel.ThreadPool'),'fbts:ProcessWorkersRequired','Use a process pool for the FDTD MEX.');
        batchTiming.worker_startup_seconds = batchTiming.worker_startup_seconds+toc(poolTimer);
        width = min(batchOptions.maxWorkers,pool.NumWorkers);
    else
        width = 1;
    end
    fbtsDirectory = fileparts(mfilename('fullpath'));
    repository = fileparts(fileparts(fileparts(fbtsDirectory)));
    workerPaths = {fbtsDirectory,fullfile(repository,'buildLib'),fullfile(repository,'forward_solver','mex')};
    tasks = cell(1,numel(pending));
    for t = 1:numel(pending)
        index = pending(t); entryIndex = find(strcmp({catalogs{j}.class_id},char(rows(index).class_id)),1);
        tasks{t} = struct('root',char(root),'row',rows(index), ...
            'entry',catalogs{j}(entryIndex),'options',options,'batchOptions',batchOptions, ...
            'references',references,'batchIdentity',identity, ...
            'codeHash',provenance.code_hash,'workerPaths',{workerPaths});
    end
    shared = struct('cfg',scenes{j},'measurements',measurements);
    runFbtsQueue(tasks,shared,width,@completed);

end

manifest.invocation_wall_seconds(end+1) = toc(timer);
persist();
summary = struct2table(rows);
fprintf('FBTS batch: %d successful, %d failed, %d pending.\n', ...
    nnz(summary.status=="success"),nnz(summary.status=="failed"),nnz(summary.status=="queued" | summary.status=="running"));

    function completed(taskIndex,row)
        rows(pending(taskIndex)) = row;
        attempted = attempted+1;
        persist();
    end

    function persist()
        writeTimer = tic;
        manifest.latest_invocation_timing = batchTiming;
        manifest.status_rows = rows;
        saveFbtsAtomic(manifestFile,struct('manifest',manifest));
        temporary = [tempname(char(root)) '.csv'];
        writetable(struct2table(rows),temporary);
        [ok,message] = movefile(temporary,fullfile(root,'cases.csv'),'f');
        assert(ok,'fbts:CommitFailed','%s',message);
        batchTiming.index_write_seconds = batchTiming.index_write_seconds+toc(writeTimer);
    end
end

function row = workerFailure(row,exception)
row.status = "failed"; row.failure_stage = "worker_or_status_write";
row.error_id = string(exception.identifier); row.error_message = string(exception.message);
end

function row = runFbtsBatchCase(task)
%runFbtsBatchCase One independent worker; writes only its own case directory.
% Process workers on the same server use the actual repository and MEX,
% rather than incomplete automatically attached dependency copies.
for k = 1:numel(task.workerPaths), addpath(task.workerPaths{k},'-begin'); end
timer = tic;
row = task.row;
directory = fullfile(task.root,char(row.output_directory));
if ~isfolder(directory), mkdir(directory); end
statusFile = fullfile(directory,'status.mat');
row.status = "running"; row.started_utc = utcNow();
saveFbtsAtomic(statusFile,struct('status',row));
stage = 'setup'; stageTimer = tic;
try
    cfg = subsetFbtsConfig(task.cfg,task.entry.active_indices);
    options = task.options;
    options.outputDirectory = "";
    execution = struct('checkpointFile',fullfile(directory,'checkpoint.mat'), ...
        'historyFile',fullfile(directory,'history_work.mat'), ...
        'checkpointEvery',task.batchOptions.checkpointEvery, ...
        'resume',task.batchOptions.resume,'computeImageError',false, ...
        'expectedCodeHash',task.codeHash);
    row.setup_seconds = toc(stageTimer);
    stage = 'fbts'; stageTimer = tic;
    [results,cfg] = runFbts(cfg,options,task.measurements,execution);
    row.fbts_seconds = toc(stageTimer);
    row.inversion_seconds = results.total_runtime;
    row.completed_iterations = results.num_iterations;
    results.timing.fbts_seconds = row.fbts_seconds;
    results.timing.setup_seconds = row.setup_seconds;
    results.timing.memory = fbtsMemorySnapshot();
    results.timing.iteration_seconds = results.iteration_runtime;
    results.timing.shared_acquisition_file = task.references.acquisition_file;
    row = solverStatistics(row,results.timing.fdtd);
    originalToLocal = zeros(1,task.cfg.antennas.numAntennas);
    originalToLocal(task.entry.active_indices) = 1:numel(task.entry.active_indices);
    caseInfo = struct('scene_id',char(row.scene_id),'class_id',char(row.class_id), ...
        'entry',task.entry,'original_to_local',originalToLocal, ...
        'tx_order_original',task.entry.active_indices,'first_transmitter',row.first_transmitter, ...
        'references',task.references,'dataset_root_relative',fullfile('..','..','..','..'), ...
        'measurement_fingerprint',task.measurements.metadata.fingerprint, ...
        'batch_identity',task.batchIdentity,'code_hash',results.code_hash);
    stage = 'save'; stageTimer = tic;
    resultFile = fullfile(task.root,char(row.result_file));
    saveFbtsBatchCase(resultFile,cfg,results,caseInfo,options);
    row.save_seconds = toc(stageTimer);
    hashTimer = tic;
    row.result_hash = string(fbtsHash(resultFile,'file'));
    row.hash_seconds = toc(hashTimer);
    row.status = "success";
    if task.batchOptions.saveFigures
        plotTimer = tic;
        try
            plotCase(directory,cfg,results);
        catch plotException
            row.plot_error = string(plotException.message);
        end
        row.plot_seconds = toc(plotTimer);
    end
catch exception
    field = [stage '_seconds'];
    if isfield(row,field) && isnan(row.(field)), row.(field) = toc(stageTimer); end
    row.status = "failed"; row.failure_stage = string(stage);
    row.error_id = string(exception.identifier); row.error_message = string(exception.message);
    checkpointFile = fullfile(directory,'checkpoint.mat');
    if isfile(checkpointFile)
        try
            cp = load(checkpointFile,'checkpoint');
            row.completed_iterations = cp.checkpoint.completed_iterations;
        catch
            % Preserve the original failure, even if a checkpoint is unreadable.
        end
    end
    failure = struct('status',row,'stack',exception.stack,'checkpoint_file','checkpoint.mat');
    if isfile(fullfile(directory,'iteration_failure.mat'))
        failure.iteration_context_file = 'iteration_failure.mat';
    end
    try
        saveFbtsAtomic(fullfile(directory,'failure.mat'),struct('failure',failure));
    catch diagnostic
        warning('fbts:FailureDiagnosticSaveFailed','%s',diagnostic.message);
    end
end
row.run_total_seconds = toc(timer); row.finished_utc = utcNow();
saveFbtsAtomic(statusFile,struct('status',row));
if row.status == "success"
    for name = {'checkpoint.mat','history_work.mat','failure.mat','iteration_failure.mat'}
        path = fullfile(directory,name{1});
        if isfile(path), delete(path); end
    end
end
end

function row = solverStatistics(row,timing)
for phase = fieldnames(timing).'
    v = timing.(phase{1}); v = v(:);
    values = [0,nan(1,7),0];
    if ~isempty(v)
        values = [numel(v),v(1),min(v),max(v),mean(v),median(v),std(v,1),var(v,1),sum(v)];
    end
    names = {'count','first_seconds','min_seconds','max_seconds','mean_seconds', ...
        'median_seconds','std_seconds','variance_seconds2','total_seconds'};
    for k = 1:numel(names), row.(['fdtd_' phase{1} '_' names{k}]) = values(k); end
end
end

function value = utcNow()
value = string(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

function plotCase(directory,cfg,results)
folder = fullfile(directory,'figures');
if ~isfolder(folder), mkdir(folder); end
for k = 1:3
    fig = figure('Visible','off'); cleanup = onCleanup(@() close(fig));
    if k < 3
        data = results.epsr_true;
        if k == 2, data = results.epsr_est; end
        imagesc(cfg.grid.originPhysical(1)+(0:cfg.Nx-1)*cfg.dx, ...
            cfg.grid.originPhysical(2)+(0:cfg.Ny-1)*cfg.dy,data.');
        axis image; set(gca,'YDir','normal'); colorbar;
        limits = [min([results.epsr_true(:);results.epsr_est(:)]), ...
            max([results.epsr_true(:);results.epsr_est(:)])];
        if limits(2)>limits(1), clim(limits); end
        xlabel('x (m)'); ylabel('y (m)');
    else
        plot(results.history.evaluated_state_index,results.cost); grid on;
        xlabel('Evaluated estimate state'); ylabel('Weighted data cost');
    end
    names = {'truth','reconstruction','cost'};
    title(names{k});
    exportgraphics(fig,fullfile(folder,[names{k} '.png']));
    clear cleanup
end
end

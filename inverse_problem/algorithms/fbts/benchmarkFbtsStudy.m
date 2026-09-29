function report = benchmarkFbtsStudy(prepared,directory,includeCuda)
%benchmarkFbtsStudy Measured, same-machine settings and 348-case estimates.
if nargin < 3, includeCuda=false; end
if ~isfolder(directory), mkdir(directory); end
provenance=fbtsCodeProvenance();
assert(strcmp(prepared.code_hash,provenance.code_hash),'fbts:StudyMismatch','Re-prepare after code/build changes.');
hardware=fbtsHardwareInfo(); timer=tic;
if includeCuda
    assert(hardware.gpu_available && hardware.cuda_mex_available,'fbts:CudaUnavailable', ...
        'CUDA benchmarking must run on the NVIDIA PC with fdtd_cuda built.');
end
report=struct('passed',false,'code_hash',provenance.code_hash,'hardware',hardware, ...
    'started_utc',fbtsUtcNow(),'screen',[],'cpu',[],'cuda',[],'recommended_options',struct());
cfg=prepared.scenes{1}; measured=scan(cfg,'cpu',1);
cfg=subsetFbtsConfig(cfg,[2 4 6 8]);
options=struct('numIterations',1,'sensitivityDownsampleFactor',prepared.sensitivity.factor);
[x,y]=find(cfg.antennas.doiMask);
% Six simultaneous DOI histories/temporaries plus process overhead; conservative.
workerBytes=6*8*(max(x)-min(x)+1)*(max(y)-min(y)+1)*cfg.Nt+1024^3;
cores=hardware.logical_processors;
if isfield(hardware,'physical_cores'), cores=hardware.physical_cores; end
budget=1;
if isfield(hardware,'MemAvailable_bytes')
    assert(workerBytes <= .85*hardware.MemAvailable_bytes,'fbts:InsufficientBenchmarkMemory', ...
        'The conservative single-worker memory estimate exceeds 85%% of available RAM.');
end
threads=2.^(0:floor(log2(min(cores,16))));
% Unmeasured warm-up, followed by bounded worker/thread configurations.
warmup=timeFbtsCase(cfg,options,measured);
if isfinite(warmup.memory.peak_resident_bytes), workerBytes=max(workerBytes,warmup.memory.peak_resident_bytes); end
if isfield(hardware,'MemAvailable_bytes'), budget=max(1,floor(.7*hardware.MemAvailable_bytes/workerBytes)); end
report.warmup=warmup; report.worker_memory_budget_bytes=workerBytes;
widths=1;
if fbtsParallelAvailable(), widths=2.^(0:floor(log2(min(cores,budget)))); end
for w=widths
    for t=threads(threads*w<=cores)
        options.solverThreads=t;
        sample=screen(w,options);
        report.screen=[report.screen;sample]; %#ok<AGROW>
        persist();
    end
end
[~,order]=sort([report.screen.cases_per_second],'descend');
for rank=1:min(2,numel(order))
    selected=report.screen(order(rank)); values=zeros(1,3); repeated=cell(1,3);
    for repetition=1:3
        options.solverThreads=selected.threads;
        sample=screen(selected.workers,options); values(repetition)=sample.cases_per_second;
        repeated{repetition}=sample;
    end
    report.screen(order(rank)).cases_per_second=median(values);
    report.screen(order(rank)).repeated_rates=values;
    report.screen(order(rank)).repeated_samples=repeated;
    persist();
end
% Select only among the two configurations actually repeated three times.
finalists=order(1:min(2,numel(order)));
rates=[report.screen(finalists).cases_per_second]; acceptable=finalists(rates>=.95*max(rates));
[~,least]=min([report.screen(acceptable).workers].*[report.screen(acceptable).threads]);
chosen=report.screen(acceptable(least));
report.recommended_options=struct('solverBackend','cpu','solverThreads',chosen.threads,'maxWorkers',chosen.workers);
report.cpu_validation=validateFbtsCpu(prepared,fullfile(directory,'cpu_validation'),chosen.threads);
report.cpu=measureBackend('cpu',chosen.threads,chosen.workers);
if includeCuda
    report.cuda_validation=validateFbtsCuda(prepared,fullfile(directory,'cuda_validation'));
    report.cuda=measureBackend('cuda',1,1);
    if report.cuda.estimated_study_seconds < .95*report.cpu.estimated_study_seconds
        report.recommended_options=struct('solverBackend','cuda','solverThreads',1,'maxWorkers',1);
    end
end
report.passed=true; report.wall_seconds=toc(timer); report.finished_utc=fbtsUtcNow(); persist();

    function measurements=scan(c,backend,scene)
        c.solverBackend=backend;
        file=fullfile(directory,sprintf('acquisition_%s_%02d.mat',backend,scene));
        if isfile(file)
            loaded=load(file,'measurements'); measurements=loaded.measurements;
            selectFbtsMeasurements(measurements,c);
        else
            measurements=acquireFbtsMeasurements(c);
            saveFbtsAtomic(file,struct('measurements',measurements));
        end
    end
    function value=screen(workers,algorithm)
        if workers>1
            pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes',max(widths)); end
            assert(~isa(pool,'parallel.ThreadPool'),'fbts:ProcessWorkersRequired','Use process workers.');
            workers=min(workers,pool.NumWorkers);
        end
        start=tic; samples=cell(1,workers);
        if workers==1
            samples{1}=timeFbtsCase(cfg,algorithm,measured);
        else
            futures=parallel.FevalFuture.empty;
            for worker=1:workers, futures(worker)=parfeval(pool,@timeFbtsCase,1,cfg,algorithm,measured); end
            cleanup=onCleanup(@() cancel(futures)); %#ok<NASGU>
            for worker=1:workers, samples{worker}=fetchOutputs(futures(worker)); end
        end
        elapsed=toc(start);
        value=struct('workers',workers,'threads',algorithm.solverThreads, ...
            'wall_seconds',elapsed,'cases_per_second',workers/elapsed,'samples',{samples}, ...
            'repeated_rates',[],'repeated_samples',{{}});
    end
    function result=measureBackend(backend,threadCount,workerCount)
        sceneNumbers=[1 4 7 10]; cases=struct([]); acquisitionSeconds=zeros(4,1);
        for shape=1:4
            c=prepared.scenes{sceneNumbers(shape)}; c.solverBackend=backend;
            measurements=scan(c,backend,sceneNumbers(shape));
            acquisitionSeconds(shape)=measurements.timing.wall_seconds;
            catalog=buildFbtsAntennaCatalog(c,'dihedral',[],123,sprintf('scene_%06d',sceneNumbers(shape)));
            half=catalog(find([catalog.disabled_count]==4,1)).active_indices;
            for active={1:8,half}
                cases(end+1)=fullCase(c,measurements,active{1},sceneNumbers(shape),backend,threadCount,numel(cases)+1); %#ok<AGROW>
            end
        end
        activeCounts=[cases.active_count]'; durations=[cases.wall_seconds]';
        coefficients=fit(activeCounts(1:6),durations(1:6));
        predicted=[ones(2,1),activeCounts(7:8)]*coefficients;
        holdout=max(abs(predicted-durations(7:8))./durations(7:8));
        if holdout>.15
            for m=[2 6]
                entry=catalog(find([catalog.disabled_count]==8-m,1));
                cases(end+1)=fullCase(c,measurements,entry.active_indices,sceneNumbers(end),backend,threadCount,numel(cases)+1); %#ok<AGROW>
            end
        end
        coefficients=fit([cases.active_count]',[cases.wall_seconds]');
        repeats=zeros(1,2); [~,middle]=sort([cases.wall_seconds]);
        for repeatCase=1:2
            entry=cases(middle(floor(numel(middle)/2)+repeatCase-1));
            c=prepared.scenes{entry.scene_index}; c.solverBackend=backend;
            measurements=scan(c,backend,entry.scene_index);
            repeatedCase=fullCase(c,measurements,entry.active_indices,entry.scene_index,backend,threadCount,numel(cases)+repeatCase);
            repeats(repeatCase)=abs(repeatedCase.wall_seconds-entry.wall_seconds)/entry.wall_seconds;
        end
        % Short single-TX check of fixed overhead, outside the fitted full runs.
        singleOptions=struct('numIterations',1,'solverBackend',backend,'solverThreads',threadCount, ...
            'sensitivityDownsampleFactor',prepared.sensitivity.factor);
        single=timeFbtsCase(subsetFbtsConfig(c,6),singleOptions,measurements);
        allDurations=[];
        for scene=1:12
            entries=buildFbtsAntennaCatalog(prepared.scenes{scene},'dihedral',[],123,sprintf('scene_%06d',scene));
            m=arrayfun(@(e) numel(e.active_indices),entries)';
            allDurations=[allDurations;max(0,[ones(numel(m),1),m]*coefficients)]; %#ok<AGROW>
        end
        % Measure concurrency with the same short workload used for tuning.
        concurrency=1; concurrent=[];
        if workerCount>1
            algorithm=singleOptions; algorithm.solverBackend='cpu'; algorithm.solverThreads=threadCount;
            serialSample=screen(1,algorithm); concurrent=screen(workerCount,algorithm);
            concurrency=concurrent.cases_per_second/serialSample.cases_per_second;
        end
        loads=zeros(1,workerCount);
        % Include per-scene barriers, just as runFbtsBatch currently does.
        scheduled=0;
        for scene=1:12
            loads(:)=0;
            for job=(scene-1)*29+(1:29)
                [~,worker]=min(loads);
                loads(worker)=loads(worker)+allDurations(job)*workerCount/concurrency;
            end
            scheduled=scheduled+max(loads);
        end
        acquisition=12*mean(acquisitionSeconds);
        uncertainty=max([.1,holdout,repeats]);
        result=struct('cases',cases,'singleton',single,'concurrent_sample',concurrent, ...
            'coefficients_seconds',coefficients,'heldout_relative_error',holdout,'repeat_relative_changes',repeats, ...
            'acquisition_seconds_per_sample',acquisitionSeconds,'estimated_serial_seconds',sum(allDurations)+acquisition, ...
            'estimated_study_seconds',scheduled+acquisition, ...
            'estimated_range_seconds',(scheduled+acquisition)*[max(0,1-uncertainty),1+uncertainty], ...
            'uncertainty_note','Observed variation/error envelope, not a statistical confidence interval.', ...
            'workers',workerCount,'threads',threadCount,'backend',backend);
    end
    function entry=fullCase(c,measurements,active,scene,backend,threadCount,number)
        algorithm=struct('numIterations',15,'solverBackend',backend,'solverThreads',threadCount, ...
            'sensitivityDownsampleFactor',prepared.sensitivity.factor);
        c=subsetFbtsConfig(c,active); startedUtc=fbtsUtcNow(); start=tic;
        file=fullfile(directory,sprintf('%s_case_%02d.mat',backend,number));
        execution=struct('computeImageError',false,'checkpointFile',[file '.checkpoint.mat'], ...
            'historyFile',[file '.history.mat'],'resume',false);
        result=runFbts(c,algorithm,measurements,execution);
        inversion=toc(start); saveTimer=tic;
        saveFbtsAtomic(file,struct('result',result,'cfg',fbtsPortableConfig(c)));
        saving=toc(saveTimer);
        hashTimer=tic; resultHash=fbtsHash(file,'file'); hashing=toc(hashTimer);
        entry=struct('scene_index',scene,'active_count',numel(active),'active_indices',active, ...
            'started_utc',startedUtc,'finished_utc',fbtsUtcNow(), ...
            'wall_seconds',inversion+saving+hashing,'inversion_seconds',inversion,'save_seconds',saving, ...
            'hash_seconds',hashing,'result_hash',resultHash, ...
            'memory',fbtsMemorySnapshot(),'result_file',file);
        saveFbtsAtomic([file '.timing.mat'],struct('entry',entry));
    end
    function persist()
        saveFbtsAtomic(fullfile(directory,'benchmark_report.mat'),struct('report',report));
    end
end

function coefficients=fit(x,y)
coefficients=[ones(numel(x),1),x]\y;
% Fixed cost and per-transmitter cost cannot be physically negative.
if any(coefficients<0)
    coefficients=[0;sum(x.*y)/sum(x.^2)];
end
end

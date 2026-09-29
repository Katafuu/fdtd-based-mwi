function report = validateFbtsCpu(prepared,directory,threads)
%validateFbtsCpu Exact three-iteration gates, including all target shapes.
if nargin < 3, threads=1; end
if ~isfolder(directory), mkdir(directory); end
provenance=fbtsCodeProvenance();
assert(strcmp(prepared.code_hash,provenance.code_hash),'fbts:StudyMismatch','Prepared code differs.');
report=struct('passed',false,'threads',threads,'code_hash',provenance.code_hash, ...
    'hardware',fbtsHardwareInfo(),'cases',[]);
ids=[1 1 1 4 7 10]; sets={1:8,[2 4 6 8],6,6,6,6};
fields={'Ez_measured','rx_model_history','raw_gradient_epsr','gradient_epsr', ...
    'direction_epsr','rx_sensitivity_coarse_history','epsr_history_doi','history'};
for k=1:numel(ids)
    cfg=subsetFbtsConfig(prepared.scenes{ids(k)},sets{k});
    options=struct('numIterations',3,'solverThreads',1, ...
        'sensitivityDownsampleFactor',prepared.sensitivity.factor);
    measured=acquireFbtsMeasurements(cfg);
    [reference,referenceStates]=runFbtsValidationCase(cfg,options,measured, ...
        fullfile(directory,sprintf('control_%02d_serial',k)));
    options.solverThreads=threads;
    [actual,actualStates]=runFbtsValidationCase(cfg,options,measured, ...
        fullfile(directory,sprintf('control_%02d_threaded',k)));
    assert(isequal(referenceStates,actualStates),'fbts:CpuMismatch','Optimizer state trajectory differs.');
    for field=fields
        assert(isequal(reference.(field{1}),actual.(field{1})), ...
            'fbts:CpuMismatch','CPU control %d differs in %s.',k,field{1});
    end
    comparison=struct('scene_index',ids(k),'active_indices',sets{k},'exact',true);
    report.cases=[report.cases;comparison]; %#ok<AGROW>
    saveFbtsAtomic(fullfile(directory,sprintf('control_%02d.mat',k)), ...
        struct('reference',reference,'actual',actual,'comparison',comparison));
end
report.passed=true;
saveFbtsAtomic(fullfile(directory,'cpu_validation.mat'),struct('report',report));
end

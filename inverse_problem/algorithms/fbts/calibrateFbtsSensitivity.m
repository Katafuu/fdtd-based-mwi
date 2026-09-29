function report = calibrateFbtsSensitivity(scenes,directory,solverThreads)
%calibrateFbtsSensitivity Freeze the largest factor passing all four shapes.
if nargin < 3, solverThreads = 1; end
if ~isfolder(directory), mkdir(directory); end
timer = tic; factors = [1 2 4]; sceneIds = [1 4 7 10];
controls=scenes(sceneIds);
for center={[200 200],[264 200]}
    builder=scenes{1}.generation.builder_options;
    builder.targetSpecs.properties.center=center{1};
    controls{end+1}=buildFbtsConfig(builder); %#ok<AGROW>
end
controlNames={'circle_random','triangle_random','square_random','hexagon_random','circle_center','circle_edge'};
provenance = fbtsCodeProvenance();
report = struct('passed',false,'factor',NaN,'factors',factors, ...
    'scene_indices',sceneIds,'controls',{controlNames},'relative_trace_error',nan(6,3), ...
    'relative_step_error',nan(6,3),'run_seconds',nan(6,3), ...
    'threshold',.05,'code_hash',provenance.code_hash);
for j = 1:numel(controls)
    cfg = controls{j}; cfg.solverThreads = solverThreads;
    scanPath = fullfile(directory,sprintf('%s_acquisition.mat',controlNames{j}));
    if isfile(scanPath)
        loaded = load(scanPath,'measurements'); measurements = loaded.measurements;
        selectFbtsMeasurements(measurements,cfg);
    else
        measurements = acquireFbtsMeasurements(cfg);
        saveFbtsAtomic(scanPath,struct('measurements',measurements));
    end
    for k = 1:3
        path = fullfile(directory,sprintf('%s_factor_%d.mat',controlNames{j},factors(k)));
        identity = fbtsHash({fbtsPortableConfig(cfg),factors(k),provenance.code_hash});
        if isfile(path)
            loaded = load(path); assert(strcmp(loaded.identity,identity),'fbts:CalibrationMismatch','Saved sensitivity calibration differs.');
            sample = loaded.sample;
        else
            runTimer = tic;
            result = runFbts(cfg,struct('numIterations',1,'sensitivityDownsampleFactor',factors(k), ...
                'solverThreads',solverThreads),measurements,struct('computeImageError',false));
            sample = struct('traces',result.Ez_sensitivity,'alpha',result.alpha, ...
                'raw_gradient',result.raw_gradient_epsr,'direction',result.direction_epsr, ...
                'step_a',result.step_a,'step_q',result.step_q,'wall_seconds',toc(runTimer));
            saveFbtsAtomic(path,struct('identity',identity,'sample',sample));
        end
        if k == 1, reference = sample; end
        report.relative_trace_error(j,k) = norm(sample.traces(:)-reference.traces(:))/max(norm(reference.traces(:)),realmin);
        report.relative_step_error(j,k) = abs(sample.alpha-reference.alpha)/max(abs(reference.alpha),realmin);
        report.run_seconds(j,k) = sample.wall_seconds;
    end
    report.wall_seconds = toc(timer);
    saveFbtsAtomic(fullfile(directory,'sensitivity_calibration.mat'),struct('report',report));
end
accepted = all(report.relative_trace_error <= .05 & report.relative_step_error <= .05,1);
report.factor = max(factors(accepted)); report.passed = ~isempty(report.factor);
assert(report.passed,'fbts:SensitivityCalibrationFailed','No sensitivity factor passed.');
report.wall_seconds = toc(timer);
saveFbtsAtomic(fullfile(directory,'sensitivity_calibration.mat'),struct('report',report));
end

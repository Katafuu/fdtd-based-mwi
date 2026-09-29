function report = calibrateFbtsStudyWindow(directory, solverThreads)
%calibrateFbtsStudyWindow Verify raw/scattered RX tails before preparing a study.
% Individual scans are saved immediately and can be reused after interruption.
if ~isfolder(directory), mkdir(directory); end
if nargin < 2, solverThreads = 1; end
timer = tic; windows = [12 24 48 96]*1e-9; scans = cell(4,3);
report = struct('passed',false,'recording_seconds',NaN,'threshold',1e-4, ...
    'significant_energy_fraction',1e-12,'windows_seconds',windows, ...
    'max_tail_fraction',nan(4,4),'max_extension_fraction',nan(4,4), ...
    'note','96 ns is diagnostic extension only; accepted recording window never exceeds 48 ns.');
for w = 1:4
    for c = 1:3
        center = [200 200]; if c == 3, center = [264 200]; end
        target = struct('name','circle','properties',struct('center',center,'radius',10), ...
            'material',struct('epsr',2,'cond_e',0,'cond_m',0,'murx',1,'mury',1));
        cfg = buildFbtsConfig(struct('seed',42,'backgroundPermittivity',45, ...
            'doiDiameter',.15,'deltaF',1/windows(w),'targetSpecs',target));
        cfg.solverThreads = solverThreads;
        if c == 1, cfg.grid.epsr = cfg.grid.epsr_bg; end
        file = fullfile(directory,sprintf('window_%dns_control_%d.mat',round(windows(w)*1e9),c));
        if isfile(file)
            saved = load(file,'measurements'); measurements = saved.measurements;
            selectFbtsMeasurements(measurements,cfg);
        else
            measurements = acquireFbtsMeasurements(cfg);
            saveFbtsAtomic(file,struct('cfg',fbtsPortableConfig(cfg),'measurements',measurements));
        end
        scans{w,c} = measurements.rx;
    end
    n = size(scans{w,1},3); tail = floor(.8*n)+1:n;
    traces = {scans{w,2},scans{w,3},scans{w,2}-scans{w,1},scans{w,3}-scans{w,1}};
    for k = 1:4, report.max_tail_fraction(w,k) = fraction(traces{k},tail); end
    if w > 1
        previousN = size(scans{w-1,1},3);
        for c=1:3
            assert(isequal(scans{w-1,c},scans{w,c}(:,:,1:previousN)), ...
                'fbts:WindowPrefixMismatch','Extending the recording window changed its earlier samples.');
        end
        for k = 1:4, report.max_extension_fraction(w-1,k) = fraction(traces{k},previousN+1:n); end
        if all(report.max_tail_fraction(w-1,:) < report.threshold) && ...
                all(report.max_extension_fraction(w-1,:) < report.threshold)
            report.passed = true; report.recording_seconds = windows(w-1);
            report.dt=cfg.dt; report.Nt=previousN;
        end
    end
    report.wall_seconds = toc(timer);
    saveFbtsAtomic(fullfile(directory,'window_calibration.mat'),struct('report',report));
    fprintf('Window %g ns: tail fractions %s\n',windows(w)*1e9,mat2str(report.max_tail_fraction(w,:),4));
    if report.passed, return; end
    if w == 3 && any(report.max_tail_fraction(w,:) >= report.threshold), break; end
end
error('fbts:RecordingWindowNotValidated', ...
    'No recording window was validated through 48 ns. Inspect window_calibration.mat before preparing the study.');
end

function value = fraction(trace,indices)
assert(all(isfinite(trace),'all'),'fbts:NonfiniteMeasurements','Nonfinite calibration trace.');
energy = sum(trace.^2,3); selected = energy > max(energy,[],'all')*1e-12;
if ~any(selected,'all'), value = 0; return; end
part = sum(trace(:,:,indices).^2,3);
value = max(part(selected)./energy(selected));
end

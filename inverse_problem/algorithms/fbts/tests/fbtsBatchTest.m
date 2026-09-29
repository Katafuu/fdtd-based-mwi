classdef fbtsBatchTest < matlab.unittest.TestCase
    %fbtsBatchTest Real solver coverage for the versioned batch dataset.
    properties
        TemporaryDirectory
    end
    methods (TestClassSetup)
        function paths(testCase)
            fbts = fileparts(fileparts(mfilename('fullpath')));
            root = fileparts(fileparts(fileparts(fbts)));
            for p = {fbts,fullfile(fbts,'tests','helpers'), ...
                    fullfile(root,'buildLib'),fullfile(root,'forward_solver','mex')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(p{1}));
            end
        end
    end
    methods (TestMethodSetup)
        function temporary(testCase)
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryDirectory = f.Folder;
        end
    end
    methods (Test)
        function dihedralCountsAndUniformOrbitMembership(testCase)
            cfg = makeFbtsTestConfig(); previous = rng;
            expected = [5 29 223]; sizes = [4 8 12];
            for j = 1:3
                cfg.antennas.numAntennas = sizes(j); cfg = buildCircularAntennaArrayIdx(cfg);
                catalog = buildFbtsAntennaCatalog(cfg,'dihedral',[],42,'scene');
                testCase.verifyNumElements(catalog,expected(j));
                testCase.verifyEqual(sum([catalog.orbit_size]),2^sizes(j)-1);
                testCase.verifyEqual(catalog,buildFbtsAntennaCatalog(cfg,'dihedral',[],42,'scene'));
                for entry = catalog.'
                    sign = 1-2*entry.reflected;
                    actual = sort(mod(sign*(entry.canonical_disabled-1)+entry.rotation,sizes(j))+1);
                    testCase.verifyEqual(entry.disabled_indices,actual);
                    testCase.verifyEqual(entry.active_indices,setdiff(1:sizes(j),actual));
                end
            end
            testCase.verifyEqual(rng,previous);
            testCase.verifyEqual(fbtsAntennaClassKey(6,[1 2 4]),fbtsAntennaClassKey(6,[1 4 6]));
            testCase.verifyNotEqual(fbtsAntennaClassKey(10,[1 2 4 7]),fbtsAntennaClassKey(10,[1 2 5 7]));
            testCase.verifyNotEqual(fbtsAntennaClassKey(6,[1 2 4],'rotation'), ...
                fbtsAntennaClassKey(6,[1 4 6],'rotation'));
            cfg.antennas.numAntennas = 8; cfg = buildCircularAntennaArrayIdx(cfg);
            c = buildFbtsAntennaCatalog(cfg,'dihedral',2);
            periodic = c(arrayfun(@(e) isequal(e.canonical_gaps,[4 4]),c));
            testCase.verifyEqual(periodic.orbit_size,4);
            testCase.verifyNumElements(buildFbtsAntennaCatalog(cfg,'rotation'),35);
            testCase.verifyNumElements(buildFbtsAntennaCatalog(cfg,'none'),255);
        end
        function catalogRejectsUnorderedOrNonuniformArrays(testCase)
            cfg = makeFbtsTestConfig(); cfg.antennas.pos([1 2],:) = cfg.antennas.pos([2 1],:);
            testCase.verifyError(@() buildFbtsAntennaCatalog(cfg),'fbts:NonuniformArray');
            testCase.verifyNumElements(buildFbtsAntennaCatalog(cfg,'none'),15);
            cfg = makeFbtsTestConfig(); cfg.dy = 2*cfg.dx;
            testCase.verifyError(@() buildFbtsAntennaCatalog(cfg),'fbts:NonuniformArray');
            cfg = makeFbtsTestConfig();
            testCase.verifyError(@() buildFbtsAntennaCatalog(cfg,'dihedral',[],0,'scene',2), ...
                'fbts:CatalogTooLarge');
        end
        function allMeasurementSubsetsMatchFreshScans(testCase)
            cfg = makeFbtsTestConfig(); full = acquireFbtsMeasurements(cfg);
            for k = 1:4
                subsets = nchoosek(1:4,k);
                for row = 1:size(subsets,1)
                    reduced = subsetFbtsConfig(cfg,subsets(row,:));
                    fresh = acquireFbtsMeasurements(reduced);
                    testCase.verifyEqual(selectFbtsMeasurements(full,reduced),fresh.rx);
                end
            end
            c = cfg; c.grid.epsr(23,25) = 3;
            testCase.verifyError(@() selectFbtsMeasurements(full,c),'fbts:MeasurementMismatch');
            c = subsetFbtsConfig(cfg,[2 4]); c.antennas.pos(1,1) = c.antennas.pos(1,1)+1;
            testCase.verifyError(@() selectFbtsMeasurements(full,c),'fbts:MeasurementMismatch');
            broken = full; broken.rx(1,1,1) = broken.rx(1,1,1)+1;
            testCase.verifyError(@() selectFbtsMeasurements(broken,cfg),'fbts:InvalidMeasurements');
        end
        function cachedInversionsMatchIndependentAcquisition(testCase)
            cfg = makeFbtsTestConfig(); full = acquireFbtsMeasurements(cfg);
            for factor = [1 4]
                for selected = {1:4,[2 4],3}
                    c = subsetFbtsConfig(cfg,selected{1});
                    options = struct('numIterations',3,'sensitivityDownsampleFactor',factor);
                    direct = runFbts(c,options); cached = runFbts(c,options,full);
                    ignored = {'timing','iteration_runtime','total_runtime','average_runtime','measurement_fingerprint'};
                    testCase.verifyEqual(rmfield(direct,ignored),rmfield(cached,ignored));
                    testCase.verifyEmpty(cached.timing.fdtd.measurement);
                end
            end
        end
        function checkpointResumeKeepsConjugateGradientTrajectory(testCase)
            cfg = makeFbtsTestConfig(); full = acquireFbtsMeasurements(cfg);
            cfg = subsetFbtsConfig(cfg,[2 4]);
            options = struct('numIterations',4);
            execution = struct('checkpointFile',fullfile(testCase.TemporaryDirectory,'checkpoint.mat'), ...
                'computeImageError',false,'onIteration',@interruptAtTwo);
            testCase.verifyError(@() runFbts(cfg,options,full,execution),'fbtsTest:Interrupted');
            stored = load(execution.checkpointFile,'checkpoint');
            testCase.verifyEqual(stored.checkpoint.completed_iterations,2);
            testCase.verifyFalse(isfield(stored.checkpoint.history,'relative_error_percent'));
            diagnostic = load(fullfile(testCase.TemporaryDirectory,'iteration_failure.mat'));
            testCase.verifyEqual(diagnostic.failureContext.phase,'iteration_callback');
            testCase.verifyEqual(diagnostic.failureContext.completed_iterations,2);
            execution.resume = true; execution.onIteration = [];
            resumed = runFbts(cfg,options,full,execution);
            direct = runFbts(cfg,options,full,struct('computeImageError',false));
            ignored = {'timing','iteration_runtime','total_runtime','average_runtime'};
            testCase.verifyEqual(rmfield(resumed,ignored),rmfield(direct,ignored));
            testCase.verifyEqual(resumed.epsr_history_doi(:,1),cfg.grid.epsr_bg(cfg.antennas.doiMask));
            testCase.verifyEqual(resumed.epsr_history_doi(:,end),resumed.epsr_est(cfg.antennas.doiMask));
            testCase.verifyEqual(resumed.history.evaluated_state_index,(0:3)');
            h = matfile(fullfile(testCase.TemporaryDirectory,'history_work.mat'),'Writable',true);
            h.rx_model_history(1,1,1,1) = h.rx_model_history(1,1,1,1)+1;
            testCase.verifyError(@() runFbts(cfg,options,full,execution),'fbts:CheckpointMismatch');
        end
        function sharedDatasetIsPortableAndMetricFree(testCase)
            cfg = makeFbtsTestConfig(); base = cfg;
            fbtsOptions = struct('numIterations',2, ...
                'outputDirectory',fullfile(testCase.TemporaryDirectory,'batch'));
            fbtsBatchOptions = struct('seed',23);
            before = findall(groot,'Type','figure');
            fbts_batch;
            testCase.verifyEqual(cfg,base);
            testCase.verifyEqual(findall(groot,'Type','figure'),before);
            testCase.verifyEqual(height(batchSummary),5);
            testCase.assertEqual(batchSummary.status,repmat("success",5,1));
            testCase.verifyEqual(batchSummary.fdtd_measurement_count,zeros(5,1));
            testCase.verifyEqual(batchSummary.fdtd_forward_count,2*batchSummary.num_active_antennas);
            acquisition = load(fullfile(batchOutputDirectory,'scenes','scene_000001','acquisition.mat'));
            testCase.verifyEqual(acquisition.measurementTiming.solver_count,4);
            for k = 1:height(batchSummary)
                data = loadFbtsBatchCase(fullfile(batchOutputDirectory,batchSummary.result_file(k)));
                testCase.verifyEqual(data.truth.epsr_true,cfg.grid.epsr);
                testCase.verifyEqual(data.cfg.antennas.originalIndices,data.result.caseInfo.entry.active_indices);
                testCase.verifyEqual(data.result.caseInfo.first_transmitter,min(data.cfg.antennas.originalIndices));
                testCase.verifyEqual(data.result.epsr_history_doi(:,end), ...
                    data.result.epsr_est(data.system.doi_mask));
                testCase.verifyFalse(isfield(data.result,'relative_error_percent'));
                testCase.verifyFalse(isfield(data.result,'Ez_measured'));
                testCase.verifyFalse(isfield(data.result.history,'relative_error_percent'));
                testCase.verifyEqual(data.result.history.evaluated_state_index,[0;1]);
                testCase.verifyEqual(size(data.result.rx_model_history,4),2);
                options = data.result.parameters.fbts_options;
                rerun = runFbts(data.cfg,options,data.measurements,struct('computeImageError',false));
                testCase.verifyEqual(rerun.epsr_est,data.result.epsr_est);
            end
            resumed = runFbtsBatch(cfg,fbtsOptions,struct('seed',23,'resume',true));
            testCase.verifyEqual(resumed,batchSummary);
            delete(fullfile(batchOutputDirectory,batchSummary.output_directory(1),'status.mat'));
            [recovered,~,savedManifest] = runFbtsBatch(cfg,fbtsOptions,struct('seed',23,'resume',true));
            testCase.verifyEqual(recovered,batchSummary);
            testCase.verifyEqual(savedManifest.catalogs,batchManifest.catalogs);
            moved = fullfile(testCase.TemporaryDirectory,'moved');
            movefile(batchOutputDirectory,moved);
            data = loadFbtsBatchCase(fullfile(moved,batchSummary.result_file(1)));
            testCase.verifyEqual(data.truth.epsr_true,cfg.grid.epsr);
            oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
            rmpath(fileparts(which('fdtd_mex')));
            library = fileparts(fileparts(which('fdtdmat.createGrid')));
            rmpath(library);
            clear fdtd_mex
            analysisOnly = loadFbtsBatchCase(fullfile(moved,batchSummary.result_file(1)));
            testCase.verifyEqual(analysisOnly.rx_selected,data.rx_selected);
        end
        function multipleScenesResumeWithoutReacquiringCompletedScene(testCase)
            a = makeFbtsTestConfig(); b = a;
            b.grid.epsr(b.grid.epsr>1) = 3; b.targets.material.epsr = 3;
            options = struct('numIterations',1,'outputDirectory',fullfile(testCase.TemporaryDirectory,'multi'));
            control = struct('disabledCounts',0,'maxCases',1);
            [first,root] = runFbtsBatch({a,b},options,control);
            testCase.verifyEqual(first.status,["success";"queued"]);
            firstHash = fbtsHash(fullfile(root,'scenes','scene_000001','acquisition.mat'),'file');
            control.resume = true; control.maxCases = Inf;
            [second,~,manifest] = runFbtsBatch({a,b},options,control);
            testCase.verifyEqual(second.status,["success";"success"]);
            testCase.verifyEqual(firstHash,fbtsHash(fullfile(root,'scenes','scene_000001','acquisition.mat'),'file'));
            testCase.verifyEqual(manifest.scene_references{1}.system_file,manifest.scene_references{2}.system_file);
            testCase.verifyNotEqual(manifest.scene_references{1}.acquisition_hash,manifest.scene_references{2}.acquisition_hash);
            options.numIterations = 2;
            testCase.verifyError(@() runFbtsBatch({a,b},options,control),'fbts:BatchMismatch');
        end
        function dryRunWritesNothingAndComputesCallBudget(testCase)
            cfg = makeFbtsTestConfig(); output = fullfile(testCase.TemporaryDirectory,'dry');
            [summary,root,manifest] = runFbtsBatch(cfg,struct('numIterations',3,'outputDirectory',output),struct('dryRun',true));
            testCase.verifyEqual(height(summary),5); testCase.verifyEqual(root,"");
            testCase.verifyFalse(isfolder(output));
            testCase.verifyEqual(manifest.resource_estimates.measurement_calls,4);
            testCase.verifyEqual(manifest.resource_estimates.fine_iteration_calls,72);
            testCase.verifyEqual(manifest.resource_estimates.coarse_iteration_calls,72);
        end
        function sparseStationaryFailuresDoNotStopOtherScenes(testCase)
            blank = makeFbtsTestConfig(); blank.grid.epsr = blank.grid.epsr_bg;
            blank.targets = blank.targets([]);
            valid = makeFbtsTestConfig();
            options = struct('numIterations',1,'outputDirectory',fullfile(testCase.TemporaryDirectory,'failures'));
            [summary,root] = runFbtsBatch({blank,valid},options,struct('disabledCounts',0));
            testCase.verifyEqual(summary.status,["failed";"success"]);
            testCase.verifyEqual(summary.failure_stage(1),"fbts");
            testCase.verifyTrue(isfile(fullfile(root,summary.output_directory(1),'failure.mat')));
            testCase.verifyTrue(isfile(fullfile(root,'cases.csv')));
        end
        function builderRetainsRandomnessAndActualMasks(testCase)
            previous = rng;
            a = buildFbtsConfig(struct('seed',52)); b = buildFbtsConfig(struct('seed',52));
            testCase.verifyEqual(a,b); testCase.verifyEqual(rng,previous);
            testCase.verifyEqual(a.generation.seed,52);
            testCase.verifyEqual(a.grid.epsr(a.targets(1).mask), ...
                repmat(a.targets(1).material.epsr,nnz(a.targets(1).mask),1));
            testCase.verifyTrue(isfield(a.targets,'physical'));
        end
        function singleAntennaSingleIterationAndChangedFiles(testCase)
            cfg = makeFbtsTestConfig(); cfg.antennas.numAntennas = 1;
            cfg.antennas.txAntennas = 1; cfg = buildCircularAntennaArrayIdx(cfg);
            cfg.source.samples = zeros(1,cfg.Nt); cfg.source.location = cfg.antennas.pos;
            options = struct('numIterations',1,'outputDirectory',fullfile(testCase.TemporaryDirectory,'single'));
            [summary,root] = runFbtsBatch(cfg,options);
            testCase.assertEqual(summary.status,"success");
            testCase.verifyEqual(summary.first_transmitter,1);
            data = loadFbtsBatchCase(fullfile(root,summary.result_file));
            testCase.verifyEqual(size(data.result.rx_model_history,4),1);
            artifact = fullfile(root,'scenes','scene_000001','acquisition.mat');
            changed = load(artifact); changed.rx_full(1,1,1) = changed.rx_full(1,1,1)+1;
            saveFbtsAtomic(artifact,changed);
            testCase.verifyError(@() loadFbtsBatchCase(fullfile(root,summary.result_file)), 'fbts:DatasetMismatch');
            testCase.verifyError(@() runFbtsBatch(cfg,options,struct('resume',true)), 'fbts:DatasetMismatch');
        end
        function directoryFailureAndCoordinatorLockAreIsolated(testCase)
            cfg = makeFbtsTestConfig();
            options = struct('numIterations',1,'outputDirectory',fullfile(testCase.TemporaryDirectory,'blocked'));
            [summary,root] = runFbtsBatch(cfg,options,struct('maxCases',1));
            blocked = fullfile(root,summary.output_directory(2));
            fid = fopen(blocked,'w'); fprintf(fid,'blocked path for test'); fclose(fid);
            continued = runFbtsBatch(cfg,options,struct('resume',true));
            testCase.verifyEqual(continued.status,["success";"failed";"success";"success";"success"]);
            lock = acquireFbtsBatchLock(char(root));
            testCase.verifyError(@() acquireFbtsBatchLock(char(root)),'fbts:BatchLocked');
            clear lock
            replacement = acquireFbtsBatchLock(char(root)); %#ok<NASGU>
        end
        function processWorkersMatchSerialCases(testCase)
            testCase.assumeTrue(fbtsParallelAvailable());
            pool = gcp('nocreate');
            created = isempty(pool);
            if created, cleanup = onCleanup(@() delete(gcp('nocreate'))); end %#ok<NASGU>
            cfg = makeFbtsTestConfig();
            options = struct('numIterations',1,'outputDirectory',fullfile(testCase.TemporaryDirectory,'serial'));
            control = struct('seed',72);
            [serial,root1] = runFbtsBatch(cfg,options,control);
            options.outputDirectory = fullfile(testCase.TemporaryDirectory,'parallel');
            control.maxWorkers = 2;
            initial=control; initial.maxCases=0;
            [planned,root2]=runFbtsBatch(cfg,options,initial);
            % Five mixed-duration cases exceed the queue width. One task has
            % an intentionally blocked directory; remaining tasks must finish.
            blocked=fullfile(root2,planned.output_directory(2));
            mkdir(fileparts(blocked)); fid=fopen(blocked,'w'); fclose(fid);
            control.resume=true;
            parallel=runFbtsBatch(cfg,options,control);
            testCase.assertEqual(parallel.status(2),"failed");
            testCase.assertEqual(parallel.status([1 3:5]),repmat("success",4,1));
            delete(blocked);
            parallel=runFbtsBatch(cfg,options,control);
            testCase.assertEqual(parallel.status,repmat("success",5,1));
            for k = 1:5
                one = load(fullfile(root1,serial.result_file(k)),'epsr_est','history');
                two = load(fullfile(root2,parallel.result_file(k)),'epsr_est','history');
                testCase.verifyEqual(one,two);
            end
        end
    end
end
function interruptAtTwo(iteration)
if iteration==2, error('fbtsTest:Interrupted','Intentional interruption after a committed checkpoint.'); end
end

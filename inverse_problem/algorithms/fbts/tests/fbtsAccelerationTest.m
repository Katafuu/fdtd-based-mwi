classdef fbtsAccelerationTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function paths(testCase)
            fbts = fileparts(fileparts(mfilename('fullpath')));
            root = fileparts(fileparts(fileparts(fbts)));
            for path = {fbts,fullfile(fbts,'tests','helpers'),fullfile(root,'buildLib'),fullfile(root,'forward_solver','mex')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(path{1}));
            end
        end
    end
    methods (Test)
        function studyGeometryIsPhysicalAndReproducible(testCase)
            before = rng;
            [scenes,index] = buildFbtsStudyScenes();
            testCase.verifyEqual(rng,before);
            testCase.verifyEqual(height(index),12);
            catalog=buildFbtsAntennaCatalog(scenes{1},'dihedral',[],123,'scene_000001');
            testCase.verifyEqual(numel(catalog),29);
            testCase.verifyEqual(sum(arrayfun(@(entry) numel(entry.active_indices),catalog)),120);
            centers = zeros(12,2);
            for k = 1:12
                cfg = scenes{k}; centers(k,:) = cfg.study.center_m;
                testCase.verifyEqual(cfg.antennas.numAntennas,8);
                testCase.verifyEqual(cfg.antennas.doiDiameter_m,.15);
                testCase.verifyEqual(cfg.grid.epsr_bg,45*ones(400));
                testCase.verifyTrue(all(ismember(cfg.grid.epsr,[2 45]),'all'));
                testCase.verifyFalse(any(cfg.targets.mask & ~cfg.antennas.doiMask,'all'));
                p = cfg.targets.physical;
                if strcmp(cfg.study.shape,'circle')
                    testCase.verifyEqual(p.radius_m,.01,'AbsTol',1e-14);
                elseif strcmp(cfg.study.shape,'square')
                    testCase.verifyEqual(diff(p.bounds_m(1:2)),.02,'AbsTol',1e-14);
                else
                    testCase.verifyEqual(vecnorm(p.vertices_m-mean(p.vertices_m,1),2,2), ...
                        .01*ones(size(p.vertices_m,1),1),'AbsTol',1e-14);
                end
            end
            for first = [1 4 7 10]
                testCase.verifyEqual(size(unique(centers(first:first+2,:),'rows'),1),3);
            end
            [again,secondIndex] = buildFbtsStudyScenes();
            testCase.verifyEqual(index,secondIndex);
            testCase.verifyEqual(again{12}.grid.epsr,scenes{12}.grid.epsr);
        end
        function regionIsExactlyAFullFieldSlice(testCase)
            cfg = makeFbtsTestConfig(); cfg.Nt = 15;
            cfg.source.samples = zeros(4,15); cfg.source.samples(2,1:3) = [1 -.5 .2];
            for tr = [false true]
                full = fdtd_mex(cfg,tr);
                for region = {[9 23 11 29],[1 1 1 48],[48 48 48 48]}
                    c = cfg; c.ezRegion = region{1}; r = region{1};
                    part = fdtd_mex(c,tr);
                    testCase.verifyEqual(part.Ez,full.Ez(r(1):r(2),r(3):r(4),:));
                    testCase.verifyEqual(part.rx_signals,full.rx_signals);
                    testCase.verifyEqual(part.ezRegion,r);
                    testCase.verifyEqual(part.Nx,cfg.Nx);
                end
            end
            cfg.ezRegion = [1 49 1 2];
            testCase.verifyError(@() fdtd_mex(cfg),'fdtd_mex:InvalidEzRegion');
            cfg.ezRegion = [3 2 1 2];
            testCase.verifyError(@() fdtd_mex(cfg),'fdtd_mex:InvalidEzRegion');
            [ez,hx,hy]=fdtd_mex(cfg); % Positional API intentionally ignores ezRegion.
            cfg=rmfield(cfg,'ezRegion'); [ezFull,hxFull,hyFull]=fdtd_mex(cfg);
            testCase.verifyEqual(ez,ezFull); testCase.verifyEqual(hx,hxFull); testCase.verifyEqual(hy,hyFull);
        end
        function threadedFieldsPreserveArithmetic(testCase)
            cfg=makeFbtsTestConfig(); cfg.Nt=32;
            cfg.source.samples=zeros(4,32); cfg.source.samples(2,1:3)=[1 -.5 .2];
            for reverse=[false true]
                cfg.solverThreads=1; reference=fdtd_mex(cfg,reverse);
                for threads=[2 4]
                    cfg.solverThreads=threads;
                    try
                        actual=fdtd_mex(cfg,reverse);
                    catch exception
                        if strcmp(exception.identifier,'fdtd_mex:OpenMPUnavailable')
                            testCase.assumeFail('Serial build: OpenMP test requires an OpenMP build.');
                        end
                        rethrow(exception)
                    end
                    testCase.verifyEqual(actual.Ez,reference.Ez);
                    testCase.verifyEqual(actual.rx_signals,reference.rx_signals);
                end
            end
        end
        function highContrastTrajectoriesMatchWithThreads(testCase)
            % Small numerical controls, not substitutes for physical window calibration.
            cfg=makeFbtsTestConfig(); cfg.antennas.numAntennas=8;
            cfg.antennas.txAntennas=1:8; cfg=buildCircularAntennaArrayIdx(cfg);
            cfg.Nt=700; cfg.deltaF=1/(cfg.Nt*cfg.dt);
            cfg.source.samples=zeros(8,cfg.Nt); cfg.source.location=cfg.antennas.pos;
            cfg.grid.epsr_bg=45*ones(cfg.Nx,cfg.Ny);
            [x,y]=ndgrid(1:cfg.Nx,1:cfg.Ny);
            shapes={'circle','triangle','square','hexagon'};
            for shape=1:4
                if shape==1
                    mask=hypot(x-23,y-25)<=2;
                elseif shape==3
                    mask=abs(x-23)<=2 & abs(y-25)<=2;
                else
                    count=3; offset=pi/2;
                    if shape==4, count=6; offset=0; end
                    angle=offset+(0:count-1)*2*pi/count;
                    mask=inpolygon(x,y,23+2*cos(angle),25+2*sin(angle));
                end
                cfg.grid.epsr=cfg.grid.epsr_bg; cfg.grid.epsr(mask)=2;
                cfg.targets=struct('name',shapes{shape},'mask',mask);
                sets={6}; if shape==1, sets={1:8,[2 4 6 8],6}; end
                for active=sets
                    c=subsetFbtsConfig(cfg,active{1}); measurements=acquireFbtsMeasurements(c);
                    options=struct('numIterations',3,'solverThreads',1,'sensitivityDownsampleFactor',1);
                    reference=runFbts(c,options,measurements,struct('computeImageError',false));
                    options.solverThreads=2;
                    try
                        actual=runFbts(c,options,measurements,struct('computeImageError',false));
                    catch exception
                        if strcmp(exception.identifier,'fdtd_mex:OpenMPUnavailable')
                            testCase.assumeFail('Serial build: requires OpenMP.');
                        end
                        rethrow(exception)
                    end
                    for field={'epsr_history_doi','rx_model_history','rx_sensitivity_coarse_history', ...
                            'raw_gradient_epsr','gradient_epsr','direction_epsr','history'}
                        testCase.verifyEqual(actual.(field{1}),reference.(field{1}));
                    end
                end
            end
        end
        function validationGateCapturesEveryState(testCase)
            cfg=makeFbtsTestConfig(); cfg.antennas.numAntennas=8;
            cfg.antennas.txAntennas=1:8; cfg=buildCircularAntennaArrayIdx(cfg);
            cfg.source.samples=zeros(8,cfg.Nt); cfg.source.location=cfg.antennas.pos;
            p=fbtsCodeProvenance();
            prepared=struct('scenes',{repmat({cfg},1,12)},'code_hash',p.code_hash, ...
                'sensitivity',struct('factor',1));
            fixture=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            report=validateFbtsCpu(prepared,fixture.Folder,1);
            testCase.verifyTrue(report.passed);
            saved=load(fullfile(fixture.Folder,'control_01_serial','iteration_states.mat'));
            testCase.verifyEqual(numel(saved.states),3);
            testCase.verifyTrue(isfield(saved.states{2},'gradEpsr'));
        end
    end
end

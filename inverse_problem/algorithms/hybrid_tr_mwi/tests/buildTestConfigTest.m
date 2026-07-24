classdef buildTestConfigTest < matlab.unittest.TestCase
    %buildTestConfigTest Verifies the compact hybrid configuration contract.

    methods (TestClassSetup)
        function addRepositoryPaths(testCase)
            hybridDir = fileparts(fileparts(mfilename('fullpath')));
            workspaceRoot = fileparts(fileparts(fileparts(hybridDir)));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'buildLib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(hybridDir, 'lib')));
        end
    end

    methods (Test, TestTags = {'Unit'})
        function producesItrCompatibleSourceAndOptions(testCase)
            hybridDir = fileparts(fileparts(mfilename('fullpath')));
            run(fullfile(hybridDir, 'build_test.m'));
            expectedPulse = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);
            inactiveAntennas = setdiff( ...
                1:cfg.antennas.numAntennas, cfg.antennas.txAntennas);
            requiredOptions = {'storeFieldHistory', 'historyStride', ...
                'numIterations', 'applyTemporalWindow', ...
                'temporalWindowTau', 'normalizeTraces'};

            testCase.verifyFalse(isfield(cfg.source, 'waveform'));
            testCase.verifySize(cfg.source.samples, ...
                [cfg.antennas.numAntennas cfg.Nt]);
            testCase.verifyEqual( ...
                cfg.source.samples(cfg.antennas.txAntennas, :), ...
                expectedPulse, AbsTol=100*eps);
            testCase.verifyEqual(cfg.source.samples(inactiveAntennas, :), ...
                zeros(numel(inactiveAntennas), cfg.Nt), AbsTol=0);
            testCase.verifyTrue(all(isfield(cfg.opts, requiredOptions)));
            testCase.verifyTrue(cfg.returnEz);
            testCase.verifyGreaterThan(cfg.snapshotStride, 0);
            testCase.verifyFalse(cfg.opts.storeFieldHistory);
            testCase.verifyEqual(cfg.opts.historyStride, 1);
            testCase.verifyEqual(cfg.opts.numIterations, 1);
            testCase.verifyTrue(cfg.opts.applyTemporalWindow);
            testCase.verifyEqual(cfg.opts.temporalWindowTau, pulse_width, ...
                AbsTol=eps(pulse_width));
            testCase.verifyTrue(cfg.opts.normalizeTraces);
        end
    end
end
classdef itrRunSamplesTest < matlab.unittest.TestCase
    %itrRunSamplesTest Tests the samples-only iterative TR source contract.

    methods (TestClassSetup)
        function addRepositoryPaths(testCase)
            hybridDir = fileparts(fileparts(mfilename('fullpath')));
            algorithmsDir = fileparts(hybridDir);
            workspaceRoot = fileparts(fileparts(fileparts(hybridDir)));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'buildLib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'forward_solver', 'mex')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(algorithmsDir, 'time_reversal', 'lib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(hybridDir, 'lib')));
        end
    end

    methods (Test, TestTags = {'Unit'})
        function rejectsWrongSampleSize(testCase)
            cfg = itrRunSamplesTest.buildConfiguration();
            cfg.source.samples = cfg.source.samples(:, 1:end-1);

            testCase.verifyError(@() itr_run(cfg), ...
                'itr_run:InvalidSourceSamples');
        end

        function rejectsNonfiniteSamples(testCase)
            cfg = itrRunSamplesTest.buildConfiguration();
            cfg.source.samples(1, 1) = NaN;

            testCase.verifyError(@() itr_run(cfg), ...
                'itr_run:InvalidSourceSamples');
        end

        function rejectsEmptySamples(testCase)
            cfg = itrRunSamplesTest.buildConfiguration();
            cfg.source.samples(:) = 0;

            testCase.verifyError(@() itr_run(cfg), ...
                'itr_run:EmptySourceSamples');
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function samplesProduceFocusResult(testCase)
            testCase.assumeEqual(exist('fdtd_mex', 'file'), 3, ...
                'The FDTD MEX file is required for the smoke test.');
            cfg = itrRunSamplesTest.buildConfiguration();

            result = itr_run(cfg);

            testCase.verifyFalse(isfield(cfg.source, 'waveform'));
            testCase.verifySize(result.focusMagFrame, [cfg.Nx cfg.Ny]);
            testCase.verifyTrue(all(isfinite(result.focusMagFrame), 'all'));
            testCase.verifyGreaterThan(result.focus_mag_peak_value, 0);
        end
    end

    methods (Static, Access = private)
        function cfg = buildConfiguration()
            cfg = struct();
            cfg.c0 = 3e8;
            cfg.mu0 = 4*pi*1e-7;
            cfg.eps0 = 8.854e-12;
            cfg.Nx = 41;
            cfg.Ny = 41;
            cfg.Nt = 220;
            cfg.sizeZ = 1;
            cfg.dx = 5e-3;
            cfg.dy = cfg.dx;
            cfg.dt = 0.99 / ...
                (cfg.c0 * sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
            cfg.snapshotStart = 0;
            cfg.snapshotStride = 1;
            cfg.returnEz = true;
            cfg.returnHx = false;
            cfg.returnHy = false;
            cfg.returnRxSignals = true;

            background = struct('epsr', 1, 'murx', 1, 'mury', 1, ...
                'cond_e', 0, 'cond_m', 0);
            cfg.grid = fdtdmat.createGrid( ...
                [cfg.Nx cfg.Ny], [cfg.dx cfg.dy], background, [0 0]);

            pmlThickness = 6;
            pmlOrder = 3;
            targetReflection = 1e-8;
            eta0 = sqrt(cfg.mu0/cfg.eps0);
            sigmaMax = -(pmlOrder + 1) * log(targetReflection) / ...
                (2 * eta0 * pmlThickness * cfg.dx);
            cfg.pml = fdtdpml.build_rectangularPML( ...
                cfg.grid, pmlThickness, pmlOrder, sigmaMax, 3);
            cfg.pml.type = 'cpml';
            cfg.pml.enabled = true;
            cfg.pml.ax = 1;
            cfg.pml.ay = 1;
            cfg.pml.az = 1;
            cfg.pml.thickness = pmlThickness;

            cfg.antennas = struct( ...
                'numAntennas', 4, ...
                'txAntennas', 1, ...
                'center', [21 21], ...
                'radius', 12, ...
                'focusPadding', 4);
            cfg = buildCircularAntennaArrayIdx(cfg);

            [gridX, gridY] = ndgrid(1:cfg.Nx, 1:cfg.Ny);
            targetMask = hypot(gridX - 21, gridY - 21) <= 2;
            cfg.grid.epsr(targetMask) = 1.1;
            cfg.grid.cond_e(targetMask) = 0.01;
            cfg.targets = struct([]);

            pulseWidth = 80e-12;
            delay = 4*pulseWidth;
            cfg.source.time = (0:cfg.Nt-1) .* cfg.dt;
            rawPulse = -((cfg.source.time - delay) ./ pulseWidth^2) .* ...
                exp(-0.5 .* ((cfg.source.time - delay) ./ pulseWidth).^2);
            normalization = max(abs(rawPulse));
            cfg.source.func = @(physicalTime) ...
                -((physicalTime - delay) ./ pulseWidth^2) .* ...
                exp(-0.5 .* ((physicalTime - delay) ./ pulseWidth).^2) ./ ...
                normalization;
            cfg.source.location = cfg.antennas.pos;
            cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);
            cfg.source.samples(cfg.antennas.txAntennas, :) = ...
                cfg.source.func(cfg.source.time);

            cfg.opts = struct( ...
                'storeFieldHistory', false, ...
                'historyStride', 1, ...
                'numIterations', 1, ...
                'applyTemporalWindow', true, ...
                'temporalWindowTau', pulseWidth, ...
                'normalizeTraces', true);
            cfg.filename = '';
        end
    end
end
classdef dbimFdtdSmokeTest < matlab.unittest.TestCase
    %dbimFdtdSmokeTest Reduced end-to-end test of the DBIM script.

    methods (TestClassSetup)
        function addRepositoryPaths(testCase)
            algorithmFolder = fileparts(fileparts(mfilename('fullpath')));
            workspaceRoot = fileparts(fileparts(fileparts(algorithmFolder)));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'buildLib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'forward_solver', 'mex')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(algorithmFolder, 'lib')));
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function oneIterationProducesBoundedResult(testCase)
            testCase.assumeEqual(exist('fdtd_mex', 'file'), 3, ...
                'The FDTD MEX file is required for the smoke test.');
            algorithmFolder = fileparts(fileparts(mfilename('fullpath')));
            cfg = dbimFdtdSmokeTest.buildConfiguration();
            outputDirectory = string(tempname);
            mkdir(outputDirectory);
            testCase.addTeardown(@() rmdir(outputDirectory, 's'));
            existingRunDirectory = fullfile(outputDirectory, 'run_0001');
            mkdir(existingRunDirectory);
            markerFile = fullfile(existingRunDirectory, 'preserve.txt');
            writelines("existing run", markerFile);
            expectedRunDirectory = fullfile(outputDirectory, 'run_0002');
            dbimOpts = struct( ...
                'frequencies', 1e9, ...
                'maxIterations', 1, ...
                'relaxation', 0.05, ...
                'residualTolerance', 1e-12, ...
                'updateTolerance', 1e-12, ...
                'epsrBounds', [1 2], ...
                'sigmaBounds', [0 0.1], ...
                'epsrScale', 1, ...
                'sigmaScale', 0.01, ...
                'operatorPrecision', "double", ...
                'assemblyBlockSize', 3, ...
                'tsvdEnergy', 0.9, ...
                'tsvdInitialRank', 4, ...
                'tsvdMaximumRank', 12, ...
                'fullSvdMaxDimension', 32, ...
                'outputDirectory', outputDirectory, ...
                'plotResolution', 72, ...
                'enablePlot', false, ...
                'verbose', false); %#ok<NASGU>

            run(fullfile(algorithmFolder, 'run_dbim.m'));

            testCase.verifyTrue(isstruct(dbimResult));
            testCase.verifySize(dbimResult.estimate.epsr, [cfg.Nx cfg.Ny]);
            testCase.verifySize(dbimResult.estimate.sigma, [cfg.Nx cfg.Ny]);
            testCase.verifyTrue(all(isfinite(dbimResult.estimate.epsr(:))));
            testCase.verifyTrue(all(isfinite(dbimResult.estimate.sigma(:))));
            testCase.verifyGreaterThanOrEqual(min(dbimResult.estimate.epsr(:)), 1);
            testCase.verifyLessThanOrEqual(max(dbimResult.estimate.epsr(:)), 2);
            testCase.verifyGreaterThanOrEqual(min(dbimResult.estimate.sigma(:)), 0);
            testCase.verifyLessThanOrEqual(max(dbimResult.estimate.sigma(:)), 0.1);
            testCase.verifyLessThanOrEqual(dbimResult.finalResidual, ...
                dbimResult.history.residual(1) + 1e-12);
            testCase.verifyEqual(size(dbimResult.pairs, 1), 6);
            testCase.verifySize(dbimResult.computationTime, [4 1]);
            testCase.verifyEqual(dbimResult.computationTimeLabels, [
                "Background/current-medium FDTD simulation"
                "FFT / frequency extraction"
                "Assembly of A and b"
                "Linear-system solve"
            ]);
            testCase.verifyTrue(all(isfinite( ...
                dbimResult.computationTime(:))));
            testCase.verifyGreaterThanOrEqual( ...
                min(dbimResult.computationTime(:)), 0);
            testCase.verifyGreaterThanOrEqual( ...
                dbimResult.measurementRuntime, 0);
            testCase.verifySize(dbimResult.iterationRuntime, [1 1]);
            testCase.verifyGreaterThanOrEqual( ...
                dbimResult.iterationRuntime, 0);
            testCase.verifyGreaterThanOrEqual(dbimResult.totalRuntime, ...
                dbimResult.measurementRuntime + ...
                sum(dbimResult.iterationRuntime));
            materialPlot = fullfile(expectedRunDirectory, ...
                '01_dbim_material_reconstruction.png');
            convergencePlot = fullfile(expectedRunDirectory, ...
                '02_dbim_convergence.png');
            resultFile = fullfile(expectedRunDirectory, 'dbim_run_data.mat');
            testCase.verifyTrue(isfile(markerFile));
            testCase.verifyTrue(isfolder(expectedRunDirectory));
            testCase.verifyEqual(string(runOutputDirectory), ...
                string(expectedRunDirectory));
            testCase.verifyTrue(isfile(materialPlot));
            testCase.verifyTrue(isfile(convergencePlot));
            testCase.verifyTrue(isfile(resultFile));

            savedData = load(resultFile);
            expectedVariables = {'cfg', 'cfg_est', 'trueCfg', 'dbimOpts', ...
                'computationTime', 'computationTimeLabels', ...
                'measurementRuntime', 'iterationRuntime', 'totalRuntime', ...
                'dbimResult'};
            testCase.verifyTrue(all(isfield(savedData, expectedVariables)));
            testCase.verifyEqual(savedData.cfg.Nx, cfg.Nx);
            testCase.verifyEqual(savedData.dbimOpts.maxIterations, 1);
            testCase.verifyEqual(savedData.computationTime, ...
                dbimResult.computationTime, AbsTol=eps);
            testCase.verifyEqual(savedData.measurementRuntime, ...
                dbimResult.measurementRuntime, AbsTol=eps);
            testCase.verifyEqual(savedData.iterationRuntime, ...
                dbimResult.iterationRuntime, AbsTol=eps);
            testCase.verifyEqual(savedData.totalRuntime, ...
                dbimResult.totalRuntime, AbsTol=eps);
            testCase.verifyEqual(dbimResult.outputFiles.figures, ...
                [string(materialPlot); string(convergencePlot)]);
            testCase.verifyEqual(dbimResult.outputFiles.matFile, ...
                string(resultFile));
        end
    end

    methods (Static, Access=private)
        function cfg = buildConfiguration()
            cfg = struct();
            cfg.c0 = 3e8;
            cfg.mu0 = 4*pi*1e-7;
            cfg.eps0 = 8.854e-12;
            cfg.Nx = 41;
            cfg.Ny = 41;
            cfg.sizeZ = 1;
            cfg.Nt = 220;
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

            pulseWidth = 80e-12;
            delay = 4*pulseWidth;
            time = (0:cfg.Nt-1) .* cfg.dt;
            rawWaveform = -((time - delay) ./ pulseWidth^2) .* ...
                exp(-0.5 .* ((time - delay) ./ pulseWidth).^2);
            sourceNormalization = max(abs(rawWaveform));
            cfg.source = struct();
            cfg.source.time = time;
            cfg.source.func = @(physicalTime) ...
                -((physicalTime - delay) ./ pulseWidth^2) .* ...
                exp(-0.5 .* ((physicalTime - delay) ./ pulseWidth).^2) ./ ...
                sourceNormalization;
            cfg.source.samples = zeros(cfg.antennas.numAntennas, cfg.Nt);
            cfg.source.samples(cfg.antennas.txAntennas, :) = ...
                cfg.source.func(cfg.source.time);
            cfg.source.location = cfg.antennas.pos;
            cfg.targets = struct([]);
        end
    end
end

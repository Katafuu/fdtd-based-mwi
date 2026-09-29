classdef fdtdMexSourcesTest < matlab.unittest.TestCase
    %FDTDMEXSOURCESTEST Validate the shared sampled multi-antenna contract.

    methods (TestClassSetup)
        function addProjectPaths(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            solverRoot = fileparts(testFolder);
            workspaceRoot = fileparts(solverRoot);
            trLibrary = fullfile(workspaceRoot, 'inverse_problem', ...
                'algorithms', 'time_reversal', 'lib');
            mexFolder = getenv('FDTD_MEX_TEST_PATH');
            if isempty(mexFolder)
                mexFolder = fullfile(solverRoot, 'mex');
            end
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                mexFolder));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'buildLib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                trLibrary));
            testCase.assertEqual(exist('fdtd_mex', 'file'), 3, ...
                'Build forward_solver/mex/fdtd_mex before running this test.');
        end
    end

    methods (Test)
        function missingNestedFieldsAreRejected(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            withoutSource = rmfield(cfg, 'source');
            withoutAntennas = rmfield(cfg, 'antennas');

            testCase.verifyError(@() fdtd_mex(withoutSource), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(withoutAntennas), ...
                'fdtd_mex:InvalidConfig');
        end

        function wrongMatrixDimensionsAreRejected(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            transposedSource = cfg;
            transposedSource.source.samples = cfg.source.samples.';
            wrongAntennaRows = cfg;
            wrongAntennaRows.antennas.pos = cfg.antennas.pos(1, :);

            testCase.verifyError(@() fdtd_mex(transposedSource), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(wrongAntennaRows), ...
                'fdtd_mex:InvalidConfig');
        end

        function invalidCoordinatesAreRejected(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            noninteger = cfg;
            noninteger.antennas.pos(1, 1) = 3.5;
            outOfBounds = cfg;
            outOfBounds.antennas.pos(2, 2) = cfg.Ny + 1;

            testCase.verifyError(@() fdtd_mex(noninteger), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(outOfBounds), ...
                'fdtd_mex:InvalidConfig');
        end

        function nonfiniteSamplesAreRejected(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            cfg.source.samples(1, 3) = NaN;

            testCase.verifyError(@() fdtd_mex(cfg), ...
                'fdtd_mex:InvalidConfig');
        end

        function receiverOutputDefaultsToEmpty(testCase)
            explicitCfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            omittedCfg = rmfield(explicitCfg, 'returnRxSignals');

            explicitResult = fdtd_mex(explicitCfg);
            omittedResult = fdtd_mex(omittedCfg);

            testCase.verifyTrue(isfield(explicitResult, 'rx_signals'));
            testCase.verifyEmpty(explicitResult.rx_signals);
            testCase.verifyEmpty(omittedResult.rx_signals);
        end

        function positionalOutputsRemainUnchanged(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            cfg.returnRxSignals = true;

            [Ez, Hx, Hy] = fdtd_mex(cfg);

            testCase.verifySize(Ez, [cfg.Nx cfg.Ny]);
            testCase.verifySize(Hx, [cfg.Nx cfg.Ny - 1]);
            testCase.verifySize(Hy, [cfg.Nx - 1 cfg.Ny]);
        end

        function receiverSignalsMatchMatlab(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 9);
            cfg.returnRxSignals = true;
            cfg.source.samples = [
                0.8 -0.2 0.1 0.4 -0.1 0.2 0 0 0
                0 0.3 -0.5 0.2 0.1 -0.2 0.1 0 0
            ];

            matlabResult = fdtdMexSourcesTest.runMatlabReference(cfg);
            mexResult = fdtd_mex(cfg);

            testCase.verifySize(mexResult.rx_signals, ...
                [cfg.antennas.numAntennas cfg.Nt]);
            testCase.verifyEqual(mexResult.rx_signals, ...
                matlabResult.rx_signals, AbsTol=1e-11);
        end

        function receiverSignalsIgnoreSnapshotStride(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 11);
            cfg.snapshotStart = 2;
            cfg.snapshotStride = 3;
            cfg.returnRxSignals = true;
            cfg.source.samples = [
                0.4 -0.1 0.2 0 0.1 0 -0.1 0 0 0 0
                0 0.2 -0.3 0.1 0 0.1 0 -0.1 0 0 0
            ];

            matlabResult = fdtdMexSourcesTest.runMatlabReference(cfg);
            mexResult = fdtd_mex(cfg);

            testCase.verifySize(mexResult.rx_signals, ...
                [cfg.antennas.numAntennas cfg.Nt]);
            testCase.verifyEqual(mexResult.rx_signals, ...
                matlabResult.rx_signals, AbsTol=1e-11);
        end

        function zeroInitialStateMatchesDefault(testCase)
            defaultCfg = fdtdMexSourcesTest.smallMexConfig(2, 8);
            defaultCfg.returnHx = true;
            defaultCfg.returnHy = true;
            seededCfg = defaultCfg;
            seededCfg.init = fdtdMexSourcesTest.zeroInitialState(seededCfg);

            defaultResult = fdtd_mex(defaultCfg);
            seededResult = fdtd_mex(seededCfg);

            testCase.verifyEqual(seededResult.Ez, defaultResult.Ez, ...
                AbsTol=1e-14);
            testCase.verifyEqual(seededResult.Hx, defaultResult.Hx, ...
                AbsTol=1e-14);
            testCase.verifyEqual(seededResult.Hy, defaultResult.Hy, ...
                AbsTol=1e-14);
        end

        function initialStateContinuesNonPmlRun(testCase)
            fullCfg = fdtdMexSourcesTest.smallMexConfig(2, 9);
            fullCfg.returnHx = true;
            fullCfg.returnHy = true;
            fullCfg.returnRxSignals = true;
            fullCfg.source.samples = [
                0.6 -0.2 0.1 0.3 -0.1 0.2 0 0.1 0
                -0.1 0.4 -0.3 0.2 0.1 -0.2 0.1 0 0
            ];
            prefixCfg = fullCfg;
            prefixCfg.Nt = 4;
            prefixCfg.source.samples = fullCfg.source.samples(:, 1:4);
            prefixCfg.returnRxSignals = false;

            fullResult = fdtd_mex(fullCfg);
            prefixResult = fdtd_mex(prefixCfg);
            continuationCfg = fullCfg;
            continuationCfg.Nt = 5;
            continuationCfg.source.samples = fullCfg.source.samples(:, 5:9);
            continuationCfg.init = struct( ...
                'Ez', prefixResult.Ez, ...
                'Hx', prefixResult.Hx, ...
                'Hy', prefixResult.Hy);
            continuationResult = fdtd_mex(continuationCfg);

            testCase.verifyEqual(continuationResult.Ez, fullResult.Ez, ...
                AbsTol=1e-11);
            testCase.verifyEqual(continuationResult.Hx, fullResult.Hx, ...
                AbsTol=1e-11);
            testCase.verifyEqual(continuationResult.Hy, fullResult.Hy, ...
                AbsTol=1e-11);
            testCase.verifyEqual(continuationResult.rx_signals, ...
                fullResult.rx_signals(:, 5:9), AbsTol=1e-11);
        end

        function invalidInitialStatesAreRejected(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            validInit = fdtdMexSourcesTest.zeroInitialState(cfg);
            nonstruct = cfg;
            nonstruct.init = 0;
            nonscalar = cfg;
            nonscalar.init = repmat(validInit, 1, 2);
            partial = cfg;
            partial.init = rmfield(validInit, 'Hy');
            wrongSize = cfg;
            wrongSize.init = validInit;
            wrongSize.init.Hx = zeros(cfg.Nx, cfg.Ny);
            complexField = cfg;
            complexField.init = validInit;
            complexField.init.Ez(1, 1) = 1i;
            logicalField = cfg;
            logicalField.init = validInit;
            logicalField.init.Hy = logical(logicalField.init.Hy);
            nonfinite = cfg;
            nonfinite.init = validInit;
            nonfinite.init.Hx(1, 1) = Inf;

            testCase.verifyError(@() fdtd_mex(nonstruct), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(nonscalar), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(partial), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(wrongSize), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(complexField), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(logicalField), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(nonfinite), ...
                'fdtd_mex:InvalidConfig');
        end

        function singleAntennaRickerMatchesMatlab(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(1, 12);
            frequency = 1.5e9;
            time = (0:cfg.Nt-1) .* cfg.dt;
            cfg.source.samples = ...
                (1 - 2 .* (pi .* (frequency .* time - 1)).^2) .* ...
                exp(-(pi .* (frequency .* time - 1)).^2);
            matlabResult = fdtdMexSourcesTest.runMatlabReference(cfg);
            mexResult = fdtd_mex(cfg);

            testCase.verifyEqual(mexResult.Ez, matlabResult.Ez, ...
                AbsTol=1e-11);
            testCase.verifyEqual(matlabResult.tx_signal, ...
                cfg.source.samples, AbsTol=0);
        end

        function multipleAntennasSuperpose(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 9);
            cfg.source.samples = [
                0.8 -0.2 0.1 0.4 -0.1 0.2 0 0 0
                0 0.3 -0.5 0.2 0.1 -0.2 0.1 0 0
            ];
            firstOnly = cfg;
            firstOnly.source.samples(2, :) = 0;
            secondOnly = cfg;
            secondOnly.source.samples(1, :) = 0;
            bothResult = fdtd_mex(cfg);
            firstResult = fdtd_mex(firstOnly);
            secondResult = fdtd_mex(secondOnly);

            testCase.verifyEqual(bothResult.Ez, ...
                firstResult.Ez + secondResult.Ez, AbsTol=1e-11);
        end

        function repeatedPositionsAccumulate(testCase)
            repeated = fdtdMexSourcesTest.smallMexConfig(2, 8);
            repeated.antennas.pos(2, :) = repeated.antennas.pos(1, :);
            repeated.source.samples = [
                0.5 -0.2 0.1 0 0.2 0 0 0
                -0.1 0.3 0.4 -0.2 0 0.1 0 0
            ];
            summed = fdtdMexSourcesTest.smallMexConfig(1, 8);
            summed.antennas.pos = repeated.antennas.pos(1, :);
            summed.source.samples = sum(repeated.source.samples, 1);
            repeatedMex = fdtd_mex(repeated);
            summedMex = fdtd_mex(summed);
            repeatedMatlab = fdtdMexSourcesTest.runMatlabReference(repeated);

            testCase.verifyEqual(repeatedMex.Ez, summedMex.Ez, ...
                AbsTol=1e-11);
            testCase.verifyEqual(repeatedMex.Ez, repeatedMatlab.Ez, ...
                AbsTol=1e-11);
        end

        function mexMatchesStandaloneFileContract(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 8);
            cfg.source.samples = [
                0.4 -0.1 0.2 0 0.1 0 0 0
                0 0.2 -0.3 0.1 0 0.1 0 0
            ];
            mexResult = fdtd_mex(cfg);
            [standaloneEz, status] = ...
                fdtdMexSourcesTest.runStandalone(cfg);

            testCase.verifyEqual(status, 0);
            testCase.verifyEqual(mexResult.Ez, standaloneEz, ...
                AbsTol=1e-12);
        end

        function standaloneRejectsTransposedSourceFile(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 8);
            [~, status] = fdtdMexSourcesTest.runStandalone( ...
                cfg, cfg.source.samples.');

            testCase.verifyNotEqual(status, 0);
        end

        function matlabHelperRequiresRowsByTime(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 7);
            [sim, grid, material, antennas, opts] = ...
                fdtdMexSourcesTest.referenceInputs(cfg);

            testCase.verifyError(@() run_fdtd_muliantenna( ...
                sim, grid, material, cfg.source.samples.', antennas, opts), ...
                'run_fdtd_muliantenna:SourceSizeMismatch');
        end

        function forwardBuildersCreateActiveCircularArrays(testCase)
            [defaultCfg, largeCfg] = ...
                fdtdMexSourcesTest.buildForwardConfigurations();
        function snapshotStacksUseCanonicalDimensions(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 4);
            cfg.snapshotStride = 1;
            cfg.returnHx = true;
            cfg.returnHy = true;
            cfg.source.samples = [0.4 -0.1 0.2 0; 0 0.3 -0.2 0.1];
            result = fdtd_mex(cfg);
            reference = fdtdMexSourcesTest.runMatlabReference(cfg);

            testCase.verifySize(result.Ez, [cfg.Nx cfg.Ny cfg.Nt]);
            testCase.verifySize(result.Hx, [cfg.Nx cfg.Ny - 1 cfg.Nt]);
            testCase.verifySize(result.Hy, [cfg.Nx - 1 cfg.Ny cfg.Nt]);
            testCase.verifyEqual(result.Ez(:, :, end), reference.Ez, AbsTol=1e-11);
            testCase.verifyEqual(result.Hx(:, :, end), ...
                reference.Hx(:, 1:cfg.Ny-1), AbsTol=1e-11);
            testCase.verifyEqual(result.Hy(:, :, end), ...
                reference.Hy(1:cfg.Nx-1, :), AbsTol=1e-11);
        end

        function numericInputsConvertToDoubleWithoutChangingOrientation(testCase)
            numericCfg = fdtdMexSourcesTest.smallMexConfig(2, 5);
            [xIndex, yIndex] = ndgrid(1:numericCfg.Nx, 1:numericCfg.Ny);
            numericCfg.grid.epsr = single(1 + 0.01*xIndex + 0.001*yIndex);
            numericCfg.grid.murx = uint16(numericCfg.grid.murx);
            numericCfg.grid.mury = single(numericCfg.grid.mury);
            numericCfg.grid.cond_e = int16(numericCfg.grid.cond_e);
            numericCfg.grid.cond_m = uint8(numericCfg.grid.cond_m);
            numericCfg.antennas.pos = uint16(numericCfg.antennas.pos);
            numericCfg.source.samples = single([ ...
                0.4 -0.1 0.2 0 0.1; 0 0.3 -0.2 0.1 0]);
            doubleCfg = numericCfg;
            doubleCfg.grid.epsr = double(numericCfg.grid.epsr);
            doubleCfg.grid.murx = double(numericCfg.grid.murx);
            doubleCfg.grid.mury = double(numericCfg.grid.mury);
            doubleCfg.grid.cond_e = double(numericCfg.grid.cond_e);
            doubleCfg.grid.cond_m = double(numericCfg.grid.cond_m);
            doubleCfg.antennas.pos = double(numericCfg.antennas.pos);
            doubleCfg.source.samples = double(numericCfg.source.samples);

            expected = fdtd_mex(doubleCfg);
            actual = fdtd_mex(numericCfg);
            testCase.verifyEqual(actual.Ez, expected.Ez, AbsTol=0);
        end

        function oldLayoutsAndUnsupportedNumericKindsAreRejected(testCase)
            cfg = fdtdMexSourcesTest.smallMexConfig(2, 5);
            oldLayout = cfg;
            oldLayout.grid.epsr = cfg.grid.epsr.';
            sparseGrid = cfg;
            sparseGrid.grid.epsr = sparse(cfg.grid.epsr);
            logicalGrid = cfg;
            logicalGrid.grid.cond_e = logical(cfg.grid.cond_e);
            complexGrid = cfg;
            complexGrid.grid.murx = cfg.grid.murx + 1i;
            pmlCfg = cfg;
            pmlCfg.pml.type = 'cpml';
            pmlCfg.pml.enabled = true;
            pmlCfg.pml.condx = zeros(cfg.Nx, cfg.Ny);
            pmlCfg.pml.condy = zeros(cfg.Nx, cfg.Ny);
            pmlCfg.pml.kx = ones(cfg.Nx, cfg.Ny);
            pmlCfg.pml.ky = ones(cfg.Nx, cfg.Ny);
            oldPmlLayout = pmlCfg;
            oldPmlLayout.pml.condx = pmlCfg.pml.condx.';

            testCase.verifyError(@() fdtd_mex(oldLayout), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(sparseGrid), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(logicalGrid), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(complexGrid), ...
                'fdtd_mex:InvalidConfig');
            testCase.verifyError(@() fdtd_mex(oldPmlLayout), ...
                'fdtd_mex:InvalidConfig');
        end


            testCase.verifyEqual(defaultCfg.antennas.numAntennas, 23);
            testCase.verifySize(defaultCfg.antennas.pos, [23 2]);
            testCase.verifySize(defaultCfg.source.samples, [23 defaultCfg.Nt]);
            testCase.verifyEqual(defaultCfg.source.samples, ...
                repmat(defaultCfg.source.samples(1, :), 23, 1), AbsTol=0);
            testCase.verifyEqual(largeCfg.antennas.numAntennas, 23);
            testCase.verifySize(largeCfg.antennas.pos, [23 2]);
            testCase.verifySize(largeCfg.source.samples, [23 largeCfg.Nt]);
            testCase.verifyEqual(largeCfg.source.samples, ...
                repmat(largeCfg.source.samples(1, :), 23, 1), AbsTol=0);
        end
    end

    methods (Static, Access = private)
        function cfg = smallMexConfig(numAntennas, numTimeSteps)
            cfg = struct();
            cfg.Nx = 9;
            cfg.Ny = 8;
            cfg.sizeZ = 1;
            cfg.Nt = numTimeSteps;
            cfg.dx = 1e-3;
            cfg.dy = 1e-3;
            cfg.dt = 1e-12;
            cfg.snapshotStart = 0;
            cfg.snapshotStride = 0;
            cfg.returnEz = true;
            cfg.returnHx = false;
            cfg.returnHy = false;
            cfg.returnRxSignals = false;

            cfg.grid = struct();
            cfg.grid.epsr = ones(cfg.Nx, cfg.Ny);
            cfg.grid.murx = ones(cfg.Nx, cfg.Ny - 1);
            cfg.grid.mury = ones(cfg.Nx - 1, cfg.Ny);
            cfg.grid.cond_e = zeros(cfg.Nx, cfg.Ny);
            cfg.grid.cond_m = zeros(cfg.Nx, cfg.Ny);

            cfg.pml = struct( ...
                'type', 'none', ...
                'enabled', false, ...
                'ax', 1, ...
                'ay', 1, ...
                'az', 1);

            basePositions = [3 3; 7 5; 4 6; 6 3];
            cfg.antennas = struct();
            cfg.antennas.numAntennas = numAntennas;
            cfg.antennas.pos = basePositions(1:numAntennas, :);

            cfg.source = struct();
            cfg.source.samples = zeros(numAntennas, numTimeSteps);
        end

        function initialState = zeroInitialState(cfg)
            initialState = struct( ...
                'Ez', zeros(cfg.Nx, cfg.Ny), ...
                'Hx', zeros(cfg.Nx, cfg.Ny - 1), ...
                'Hy', zeros(cfg.Nx - 1, cfg.Ny));
        end

        function result = runMatlabReference(cfg)
            [sim, grid, material, antennas, opts] = ...
                fdtdMexSourcesTest.referenceInputs(cfg);
            result = run_fdtd_muliantenna(sim, grid, material, ...
                cfg.source.samples, antennas, opts);
        end

        function [sim, grid, material, antennas, opts] = referenceInputs(cfg)
            sim = struct( ...
                'c0', 1 / sqrt((4*pi*1e-7) * 8.854187817e-12), ...
                'mu0', 4*pi*1e-7, ...
                'eps0', 8.854187817e-12, ...
                'Nx', cfg.Nx, ...
                'Ny', cfg.Ny, ...
                'dx', cfg.dx, ...
                'dy', cfg.dy, ...
                'dt', cfg.dt, ...
                'Nt', cfg.Nt);
            grid = struct( ...
                'sizeXY', [cfg.Nx cfg.Ny], ...
                'spacingXY', [cfg.dx cfg.dy]);
            material = struct( ...
                'eps_r', cfg.grid.epsr, ...
                'sigma', cfg.grid.cond_e, ...
                'sigma_e', zeros(cfg.Nx, cfg.Ny));
            antennas = cfg.antennas;
            opts = struct('storeFieldHistory', false);
        end

        function [ez, status] = runStandalone(cfg, sourceSamples)
            if nargin < 2
                sourceSamples = cfg.source.samples;
            end

            testFolder = fileparts(mfilename('fullpath'));
            solverRoot = fileparts(testFolder);
            executable = fullfile(solverRoot, 'tmzdemo2.exe');
            inputFolder = tempname;
            mkdir(inputFolder);
            cleanup = onCleanup(@() rmdir(inputFolder, 's'));
            outputPath = fullfile(inputFolder, 'ez_field.csv');

            configLines = [
                "Nx " + cfg.Nx
                "Ny " + cfg.Ny
                "Nt " + cfg.Nt
                "dx " + compose('%.17g', cfg.dx)
                "dy " + compose('%.17g', cfg.dy)
                "dt " + compose('%.17g', cfg.dt)
                "num_antennas " + cfg.antennas.numAntennas
                "pml_type none"
                "pml_enabled 0"
                "snapshot_start " + (cfg.Nt - 1)
                "snapshot_stride 0"
            ];
            writelines(configLines, fullfile(inputFolder, 'grid_config.txt'));
            writematrix(cfg.antennas.pos, ...
                fullfile(inputFolder, 'antennas.csv'));
            writematrix(sourceSamples, ...
                fullfile(inputFolder, 'source.csv'));
            writematrix(cfg.grid.epsr, fullfile(inputFolder, 'epsr.csv'));
            writematrix(cfg.grid.murx, fullfile(inputFolder, 'murx.csv'));
            writematrix(cfg.grid.mury, fullfile(inputFolder, 'mury.csv'));
            writematrix(cfg.grid.cond_e, ...
                fullfile(inputFolder, 'cond_e.csv'));
            writematrix(cfg.grid.cond_m, ...
                fullfile(inputFolder, 'cond_m.csv'));

            command = sprintf('"%s" "%s" "%s"', ...
                executable, inputFolder, outputPath);
            [status, ~] = system(command);
            if status ~= 0
                ez = [];
                return;
            end
            snapshots = readmatrix(outputPath, 'NumHeaderLines', 1);
            ez = reshape(snapshots(end, :), [cfg.Ny cfg.Nx]).';
        end

        function [defaultCfg, largeCfg] = buildForwardConfigurations()
            testFolder = fileparts(mfilename('fullpath'));
            solverRoot = fileparts(testFolder);
            run(fullfile(solverRoot, 'examples', 'build_fdtd_case.m'));
            defaultCfg = cfg;
            clear cfg grid pml
            run(fullfile(solverRoot, 'examples', ...
                'build_large_rectangular_pml_case.m'));
            largeCfg = cfg;
        end
    end
end

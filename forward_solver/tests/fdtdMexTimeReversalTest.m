classdef fdtdMexTimeReversalTest < matlab.unittest.TestCase
    %FDTDMEXTIMEREVERSALTEST Validate loss-uncompensated MEX TR updates.

    methods (TestClassSetup)
        function addMexPath(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            solverRoot = fileparts(testFolder);
            mexFolder = getenv('FDTD_MEX_TEST_PATH');
            if isempty(mexFolder)
                mexFolder = fullfile(solverRoot, 'mex');
            end
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                mexFolder));
            testCase.assertEqual(exist('fdtd_mex', 'file'), 3, ...
                'Build forward_solver/mex/fdtd_mex before running this test.');
        end
    end

    methods (Test)
        function omittedFlagMatchesExplicitFalse(testCase)
            cfg = fdtdMexTimeReversalTest.smallConfig(2, 6);
            cfg.source.samples = [
                0.4 -0.2 0.1 0 0.05 0
                0 0.3 -0.1 0.2 0 -0.05
            ];

            omitted = fdtd_mex(cfg);
            explicit = fdtd_mex(cfg, false);

            testCase.verifyEqual(explicit.Ez, omitted.Ez, AbsTol=0);
            testCase.verifyEqual(explicit.Hx, omitted.Hx, AbsTol=0);
            testCase.verifyEqual(explicit.Hy, omitted.Hy, AbsTol=0);
            testCase.verifyEqual(explicit.rx_signals, omitted.rx_signals, ...
                AbsTol=0);
        end

        function invalidFlagsAreRejected(testCase)
            cfg = fdtdMexTimeReversalTest.smallConfig(1, 2);

            testCase.verifyError(@() fdtd_mex(cfg, 1), ...
                'fdtd_mex:InvalidInput');
            testCase.verifyError(@() fdtd_mex(cfg, logical([1 0])), ...
                'fdtd_mex:InvalidInput');
            testCase.verifyError(@() fdtd_mex(cfg, true, false), ...
                'fdtd_mex:InvalidInput');
        end

        function oneStepMatchesSignReversedReference(testCase)
            cfg = fdtdMexTimeReversalTest.smallConfig(1, 1);
            cfg.init = fdtdMexTimeReversalTest.seedInitialState(cfg);

            expected = fdtdMexTimeReversalTest.runTrReference(cfg);
            actual = fdtd_mex(cfg, true);

            fdtdMexTimeReversalTest.verifyFields(testCase, actual, expected);
        end

        function conductivityRemainsLossy(testCase)
            cfg = fdtdMexTimeReversalTest.smallConfig(1, 7);
            cfg.grid.epsr = 2.4 .* ones(cfg.Nx, cfg.Ny);
            cfg.grid.cond_e = 0.18 .* ones(cfg.Nx, cfg.Ny);
            cfg.init = fdtdMexTimeReversalTest.seedInitialState(cfg);

            expected = fdtdMexTimeReversalTest.runTrReference(cfg);
            actual = fdtd_mex(cfg, true);

            fdtdMexTimeReversalTest.verifyFields(testCase, actual, expected);
        end

        function sourceInjectionAndReceiverSamplingMatch(testCase)
            cfg = fdtdMexTimeReversalTest.smallConfig(2, 5);
            cfg.source.samples = [
                0.7 -0.2 0.1 0 0.05
                -0.1 0.4 -0.3 0.2 0
            ];

            expected = fdtdMexTimeReversalTest.runTrReference(cfg);
            actual = fdtd_mex(cfg, true);

            fdtdMexTimeReversalTest.verifyFields(testCase, actual, expected);
            testCase.verifyEqual(actual.rx_signals, expected.rx_signals, ...
                AbsTol=1e-11);
        end

        function cpmlKeepsConvolutionAndReversesCoupling(testCase)
            cfg = fdtdMexTimeReversalTest.smallConfig(1, 1);
            cfg.init = fdtdMexTimeReversalTest.seedInitialState(cfg);
            cfg.pml.type = 'cpml';
            cfg.pml.enabled = true;
            cfg.pml.ax = 0;
            cfg.pml.ay = 0;
            cfg.pml.condx = 0.025 .* ones(cfg.Nx, cfg.Ny);
            cfg.pml.condy = 0.035 .* ones(cfg.Nx, cfg.Ny);
            cfg.pml.kx = 1.2 .* ones(cfg.Nx, cfg.Ny);
            cfg.pml.ky = 1.3 .* ones(cfg.Nx, cfg.Ny);

            expected = fdtdMexTimeReversalTest.runTrReference(cfg);
            actual = fdtd_mex(cfg, true);

            fdtdMexTimeReversalTest.verifyFields(testCase, actual, expected);
        end
    end

    methods (Static, Access = private)
        function cfg = smallConfig(numAntennas, numSteps)
            cfg = struct();
            cfg.Nx = 8;
            cfg.Ny = 7;
            cfg.sizeZ = 1;
            cfg.Nt = numSteps;
            cfg.dx = 1e-3;
            cfg.dy = 1e-3;
            cfg.dt = 1e-12;
            cfg.snapshotStart = 0;
            cfg.snapshotStride = 0;
            cfg.returnEz = true;
            cfg.returnHx = true;
            cfg.returnHy = true;
            cfg.returnRxSignals = true;

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

            positions = [3 3; 6 5; 4 6];
            cfg.antennas = struct( ...
                'numAntennas', numAntennas, ...
                'pos', positions(1:numAntennas, :));
            cfg.source = struct( ...
                'samples', zeros(numAntennas, numSteps));
        end

        function init = seedInitialState(cfg)
            [xE, yE] = ndgrid(0:cfg.Nx-1, 0:cfg.Ny-1);
            [xHx, yHx] = ndgrid(0:cfg.Nx-1, 0:cfg.Ny-2);
            [xHy, yHy] = ndgrid(0:cfg.Nx-2, 0:cfg.Ny-1);

            init = struct();
            init.Ez = 0.03 .* sin(0.4 .* xE) .* cos(0.3 .* yE);
            init.Hx = 2e-5 .* cos(0.2 .* xHx + 0.3 .* yHx);
            init.Hy = 3e-5 .* sin(0.3 .* xHy + 0.2 .* yHy);
        end

        function result = runTrReference(cfg)
            eps0 = 8.854187817e-12;
            mu0 = 4*pi*1e-7;
            Nx = cfg.Nx;
            Ny = cfg.Ny;
            dt = cfg.dt;

            epsr = cfg.grid.epsr;
            murx = cfg.grid.murx;
            mury = cfg.grid.mury;
            condE = cfg.grid.cond_e;
            condM = cfg.grid.cond_m;

            epsMat = eps0 .* epsr;
            ae = condE .* dt ./ (2 .* epsMat);
            ceze = (1 - ae) ./ (1 + ae);
            cezh = (dt ./ (epsMat .* cfg.dx)) ./ (1 + ae);

            muX = mu0 .* murx;
            amX = condM(:, 1:Ny-1) .* dt ./ (2 .* muX);
            chxh = (1 - amX) ./ (1 + amX);
            chxe = (dt ./ (muX .* cfg.dy)) ./ (1 + amX);

            muY = mu0 .* mury;
            amY = condM(1:Nx-1, :) .* dt ./ (2 .* muY);
            chyh = (1 - amY) ./ (1 + amY);
            chye = (dt ./ (muY .* cfg.dx)) ./ (1 + amY);

            if isfield(cfg, 'init')
                Ez = cfg.init.Ez;
                Hx = cfg.init.Hx;
                Hy = cfg.init.Hy;
            else
                Ez = zeros(Nx, Ny);
                Hx = zeros(Nx, Ny - 1);
                Hy = zeros(Nx - 1, Ny);
            end

            [pmlEnabled, kx, ky, bx, by, cx, cy] = ...
                fdtdMexTimeReversalTest.pmlCoefficients(cfg, eps0);
            psiEzX = zeros(Nx, Ny);
            psiEzY = zeros(Nx, Ny);
            psiHxY = zeros(Nx, Ny - 1);
            psiHyX = zeros(Nx - 1, Ny);
            receiverIndex = sub2ind([Nx Ny], ...
                cfg.antennas.pos(:, 1), cfg.antennas.pos(:, 2));
            rxSignals = zeros(cfg.antennas.numAntennas, cfg.Nt);
            ii = 2:Nx-1;
            jj = 2:Ny-1;

            for n = 1:cfg.Nt
                dEzY = Ez(:, 2:Ny) - Ez(:, 1:Ny-1);
                dEzX = Ez(2:Nx, :) - Ez(1:Nx-1, :);
                Hx = chxh .* Hx + chxe .* dEzY ./ ky(:, 1:Ny-1);
                Hy = chyh .* Hy - chye .* dEzX ./ kx(1:Nx-1, :);

                if pmlEnabled
                    psiHxY = cy(:, 1:Ny-1) .* (dEzY ./ cfg.dy) + ...
                        by(:, 1:Ny-1) .* psiHxY;
                    psiHyX = cx(1:Nx-1, :) .* (dEzX ./ cfg.dx) + ...
                        bx(1:Nx-1, :) .* psiHyX;
                    Hx = Hx + (dt ./ muX) .* psiHxY;
                    Hy = Hy - (dt ./ muY) .* psiHyX;
                end

                curlH = ...
                    (Hy(ii, jj) - Hy(ii-1, jj)) ./ kx(ii, jj) - ...
                    (Hx(ii, jj) - Hx(ii, jj-1)) ./ ky(ii, jj);
                Ez(ii, jj) = ceze(ii, jj) .* Ez(ii, jj) - ...
                    cezh(ii, jj) .* curlH;

                if pmlEnabled
                    dHyX = (Hy(ii, jj) - Hy(ii-1, jj)) ./ cfg.dx;
                    dHxY = (Hx(ii, jj) - Hx(ii, jj-1)) ./ cfg.dy;
                    psiEzX(ii, jj) = cx(ii, jj) .* dHyX + ...
                        bx(ii, jj) .* psiEzX(ii, jj);
                    psiEzY(ii, jj) = cy(ii, jj) .* dHxY + ...
                        by(ii, jj) .* psiEzY(ii, jj);
                    Ez(ii, jj) = Ez(ii, jj) - ...
                        (dt ./ epsMat(ii, jj)) .* ...
                        (psiEzX(ii, jj) - psiEzY(ii, jj));
                end

                for antenna = 1:cfg.antennas.numAntennas
                    pos = cfg.antennas.pos(antenna, :);
                    Ez(pos(1), pos(2)) = Ez(pos(1), pos(2)) + ...
                        cfg.source.samples(antenna, n);
                end
                rxSignals(:, n) = Ez(receiverIndex);
            end

            result = struct('Ez', Ez, 'Hx', Hx, 'Hy', Hy, ...
                'rx_signals', rxSignals);
        end

        function [enabled, kx, ky, bx, by, cx, cy] = ...
                pmlCoefficients(cfg, eps0)
            enabled = cfg.pml.enabled && ...
                ~any(strcmp(cfg.pml.type, {'none', 'off'}));
            kx = ones(cfg.Nx, cfg.Ny);
            ky = ones(cfg.Nx, cfg.Ny);
            condx = zeros(cfg.Nx, cfg.Ny);
            condy = zeros(cfg.Nx, cfg.Ny);

            if enabled
                kx = cfg.pml.kx;
                ky = cfg.pml.ky;
                condx = cfg.pml.condx;
                condy = cfg.pml.condy;
            end

            bx = exp(-((cfg.pml.ax ./ eps0) + ...
                condx ./ (kx .* eps0)) .* cfg.dt);
            by = exp(-((cfg.pml.ay ./ eps0) + ...
                condy ./ (ky .* eps0)) .* cfg.dt);
            denomX = condx .* kx + kx.^2 .* cfg.pml.ax;
            denomY = condy .* ky + ky.^2 .* cfg.pml.ay;
            cx = zeros(cfg.Nx, cfg.Ny);
            cy = zeros(cfg.Nx, cfg.Ny);
            maskX = denomX ~= 0;
            maskY = denomY ~= 0;
            cx(maskX) = condx(maskX) ./ denomX(maskX) .* (bx(maskX) - 1);
            cy(maskY) = condy(maskY) ./ denomY(maskY) .* (by(maskY) - 1);
        end

        function verifyFields(testCase, actual, expected)
            testCase.verifyEqual(actual.Ez, expected.Ez, AbsTol=1e-11);
            testCase.verifyEqual(actual.Hx, expected.Hx, AbsTol=1e-11);
            testCase.verifyEqual(actual.Hy, expected.Hy, AbsTol=1e-11);
        end
    end
end

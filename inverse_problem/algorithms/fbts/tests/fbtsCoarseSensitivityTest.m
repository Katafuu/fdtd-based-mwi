classdef fbtsCoarseSensitivityTest < matlab.unittest.TestCase
    %fbtsCoarseSensitivityTest Smoke-test the coarse two-solve derivative.

    methods (TestClassSetup)
        function addRepositoryPaths(testCase)
            testFile = mfilename('fullpath');
            fbtsDir = fileparts(fileparts(testFile));
            algorithmsDir = fileparts(fbtsDir);
            repoRoot = fileparts(fileparts(algorithmsDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(repoRoot, 'buildLib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(repoRoot, 'forward_solver', 'mex')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fbtsDir));
        end
    end

    methods (Test)
        function coarseFiniteDifferenceHasConsistentTimeGrid(testCase)
            cfg = makeSmallConfig();
            perturbationStep = 1e-3;
            perturbedCfg = cfg;
            perturbedCfg.grid.epsr(7:10, 7:10) = ...
                perturbedCfg.grid.epsr(7:10, 7:10) + perturbationStep;

            coarseBaseCfg = prepareCoarseSensitivityCfg(cfg, 4);
            coarsePerturbedCfg = prepareCoarseSensitivityCfg(perturbedCfg, 4);
            coarseTime = coarseBaseCfg.source.time;
            sourcePulse = zeros(1, coarseBaseCfg.Nt);
            sourcePulse(1) = 1;
            coarseBaseCfg.source.samples(1, :) = sourcePulse;
            coarsePerturbedCfg.source.samples(1, :) = sourcePulse;

            baselineResult = fdtd_mex(coarseBaseCfg);
            perturbedResult = fdtd_mex(coarsePerturbedCfg);
            coarseSensitivity = (perturbedResult.rx_signals - ...
                baselineResult.rx_signals) ./ perturbationStep;
            fineTime = (0:cfg.Nt-1) .* cfg.dt;
            fineSensitivity = interp1( ...
                coarseTime, coarseSensitivity.', fineTime, 'linear', 0).';
            coarseWeight = ones(1, coarseBaseCfg.Nt);
            stepA = sum(coarseWeight .* coarseSensitivity.^2, 'all') .* ...
                coarseBaseCfg.dt;

            testCase.verifyEqual([coarseBaseCfg.Nx coarseBaseCfg.Ny], [4 4]);
            testCase.verifyEqual(coarseBaseCfg.dt, 4 * cfg.dt, RelTol=1e-12);
            testCase.verifyEqual(coarseBaseCfg.Nt, ...
                ceil(1 / (coarseBaseCfg.dt * cfg.deltaF)));
            testCase.verifyEqual(size(coarseBaseCfg.source.samples), ...
                [cfg.antennas.numAntennas coarseBaseCfg.Nt]);
            testCase.verifyEqual(size(coarseSensitivity), ...
                [cfg.antennas.numAntennas coarseBaseCfg.Nt]);
            testCase.verifyEqual(size(fineSensitivity), ...
                [cfg.antennas.numAntennas cfg.Nt]);
            testCase.verifyTrue(all(isfinite(coarseSensitivity), 'all'));
            testCase.verifyTrue(isfinite(stepA) && stepA >= 0);
        end
    end
end

function cfg = makeSmallConfig()
cfg = struct();
cfg.c0 = 3e8;
cfg.mu0 = 4*pi*1e-7;
cfg.eps0 = 8.854e-12;
cfg.Nx = 16;
cfg.Ny = 16;
cfg.sizeZ = 1;
cfg.dx = 1e-2;
cfg.dy = 1e-2;
cfg.dt = 1 / (cfg.c0 * sqrt(1/cfg.dx^2 + 1/cfg.dy^2));
cfg.deltaF = 1 / (32 * cfg.dt);
cfg.Nt = ceil(1 / (cfg.dt * cfg.deltaF));
cfg.snapshotStart = 0;
cfg.snapshotStride = 0;
cfg.returnEz = false;
cfg.returnHx = false;
cfg.returnHy = false;
cfg.returnRxSignals = true;
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy]);
cfg.pml = fdtdpml.build_rectangularPML(cfg.grid, 4, 2, 2, 2);
cfg.pml.type = 'cpml';
cfg.pml.enabled = true;
cfg.pml.ax = 1;
cfg.pml.ay = 1;
cfg.pml.az = 1;
cfg.pml.thickness = 4;
cfg.antennas = struct( ...
    'numAntennas', 2, ...
    'pos', [5 5; 9 9], ...
    'center', [7 7], ...
    'radius', 4, ...
    'pmlPadding', 0, ...
    'focusPadding', 0, ...
    'doiMask', true(cfg.Nx, cfg.Ny));
cfg.source = struct( ...
    'location', cfg.antennas.pos, ...
    'time', (0:cfg.Nt-1) .* cfg.dt, ...
    'samples', zeros(cfg.antennas.numAntennas, cfg.Nt));
end

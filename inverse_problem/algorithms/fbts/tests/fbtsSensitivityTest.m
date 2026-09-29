classdef fbtsSensitivityTest < matlab.unittest.TestCase
    %fbtsSensitivityTest Verify the forward signals serve as the full-grid baseline.

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
        function forwardSignalsMatchIndependentBaseline(testCase)
            cfg = makeSmallConfig();
            sourcePulse = zeros(1, cfg.Nt);
            sourcePulse(1) = 1;
            cfg.source.samples(1, :) = sourcePulse;
            cfg.returnEz = true;

            forwardResult = fdtd_mex(cfg);
            baselineCfg = cfg;
            baselineCfg.returnEz = false;
            baselineResult = fdtd_mex(baselineCfg);

            perturbationStep = 1e-3;
            perturbedCfg = baselineCfg;
            perturbedCfg.grid.epsr(7:10, 7:10) = ...
                perturbedCfg.grid.epsr(7:10, 7:10) + perturbationStep;
            perturbedResult = fdtd_mex(perturbedCfg);

            sensitivityEz = (perturbedResult.rx_signals - ...
                forwardResult.rx_signals) ./ perturbationStep;
            referenceSensitivity = (perturbedResult.rx_signals - ...
                baselineResult.rx_signals) ./ perturbationStep;

            testCase.verifyEqual(forwardResult.rx_signals, ...
                baselineResult.rx_signals);
            testCase.verifyEqual(sensitivityEz, referenceSensitivity);
            testCase.verifyEqual(size(sensitivityEz), ...
                [cfg.antennas.numAntennas cfg.Nt]);
            testCase.verifyTrue(all(isfinite(sensitivityEz), 'all'));
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

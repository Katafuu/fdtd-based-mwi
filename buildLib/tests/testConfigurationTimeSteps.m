classdef testConfigurationTimeSteps < matlab.unittest.TestCase
    %testConfigurationTimeSteps Verify build scripts derive Nt from dt and deltaF.

    methods (Test)
        function inverseProblemBuildersUseMinimumTimeStepFormula(testCase)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            builderPaths = {
                fullfile(repoRoot, 'inverse_problem', 'algorithms', ...
                    'dbim', 'build_cfg_dbim.m')
                fullfile(repoRoot, 'inverse_problem', 'algorithms', ...
                    'fbts', 'build_cfg.m')
                fullfile(repoRoot, 'inverse_problem', 'algorithms', ...
                    'hybrid_tr_mwi', 'build_cfg.m')
                fullfile(repoRoot, 'inverse_problem', 'algorithms', ...
                    'hybrid_tr_mwi', 'build_test.m')
                fullfile(repoRoot, 'inverse_problem', 'algorithms', ...
                    'tds_verification', 'build_cfg.m')
                fullfile(repoRoot, 'inverse_problem', 'algorithms', ...
                    'time_reversal', 'build_test.m')
            };
            expectedNt = [4243 2547 600 800 600 800];

            for builderIndex = 1:numel(builderPaths)
                run(builderPaths{builderIndex});

                expectedTimeSteps = ceil(1 / (cfg.dt * cfg.deltaF));
                if builderIndex == 2
                    expectedTimeSteps = max(expectedTimeSteps, ...
                        ceil(6e-9 / cfg.dt) + 1);
                    testCase.verifyGreaterThanOrEqual( ...
                        (cfg.Nt - 1) * cfg.dt, 6e-9);
                end
                testCase.verifyEqual(cfg.Nt, expectedTimeSteps, ...
                    sprintf('Formula mismatch in %s.', builderPaths{builderIndex}));
                testCase.verifyEqual(cfg.Nt, expectedNt(builderIndex), ...
                    sprintf('Unexpected Nt in %s.', builderPaths{builderIndex}));
                testCase.verifyEqual(size(cfg.source.samples), ...
                    [cfg.antennas.numAntennas cfg.Nt], ...
                    sprintf('Source size mismatch in %s.', builderPaths{builderIndex}));
            end
        end
    end
end

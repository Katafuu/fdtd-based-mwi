classdef hybridMaximumChordTest < matlab.unittest.TestCase
    %hybridMaximumChordTest Confirm supplied support controls the inversion.

    methods (TestClassSetup)
        function addHybridLibrary(testCase)
            libraryDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                libraryDir));
        end
    end

    methods (Test)
        function testSelectedChordsFollowSuppliedMask(testCase)
            horizontalMask = false(32, 22);
            horizontalMask(6:26, 10:12) = true;
            verticalMask = false(32, 22);
            verticalMask(15:17, 4:18) = true;

            measurement = struct( ...
                'gridSize', [32 22], ...
                'antennaPositions', [2 11; 30 11; 16 2; 16 20], ...
                'dx', 1e-3, 'dy', 1e-3, ...
                'pairEligible', logical([0 1 0 0; 1 0 0 0; ...
                    0 0 0 1; 0 0 1 0]), ...
                'matchedDelaySteps', 10*ones(4), ...
                'c0', 3e8, 'dt', 1e-12, 'backgroundEpsr', 1);
            measurement.targetMask = verticalMask; % Must be ignored.

            [~, horizontal] = estimateMaximumChordEpsr( ...
                measurement, horizontalMask);
            [~, vertical] = estimateMaximumChordEpsr( ...
                measurement, verticalMask);

            testCase.verifyTrue(horizontal.selectedPairs(1, 2));
            testCase.verifyFalse(horizontal.selectedPairs(3, 4));
            testCase.verifyTrue(vertical.selectedPairs(3, 4));
            testCase.verifyFalse(vertical.selectedPairs(1, 2));
        end

        function testReconstructionUsesOnlySuppliedSupport(testCase)
            supportMask = false(9, 7);
            supportMask(3:5, 2:4) = true;
            reconstruction = reconstructHomogeneousTargetMap( ...
                supportMask, 1.25, 4.5);

            testCase.verifyEqual(reconstruction(supportMask), ...
                4.5*ones(nnz(supportMask), 1));
            testCase.verifyEqual(reconstruction(~supportMask), ...
                1.25*ones(nnz(~supportMask), 1));
        end
    end
end

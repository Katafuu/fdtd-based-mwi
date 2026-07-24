classdef tdVerificationHelpersTest < matlab.unittest.TestCase
    %tdVerificationHelpersTest Focused tests for the verification workflow.

    methods (TestClassSetup)
        function addVerificationLibrary(testCase)
            testDirectory = fileparts(mfilename('fullpath'));
            libraryDirectory = fullfile(fileparts(testDirectory), 'lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                libraryDirectory));
        end
    end

    methods (Test)
        function testHorizontalRayHasFifteenPixelInteriorWidth(testCase)
            rayMask = buildThickRayMask( ...
                [41 41], [10 21], [32 21], 15);

            testCase.verifyEqual(nnz(rayMask(21, :)), 15);
            testCase.verifyTrue(all(rayMask(21, 14:28)));
            testCase.verifyFalse(rayMask(21, 13));
            testCase.verifyFalse(rayMask(21, 29));
            testCase.verifyFalse(rayMask(1, 21));
            testCase.verifyFalse(rayMask(41, 21));
        end

        function testArithmeticMixtureRecoversKnownPermittivity(testCase)
            backgroundEpsr = 1.0;
            knownTargetEpsr = 5.0;
            numTargetPixels = 30;
            numBackgroundPixels = 70;
            averageEpsr = (numTargetPixels*knownTargetEpsr + ...
                numBackgroundPixels*backgroundEpsr) / ...
                (numTargetPixels + numBackgroundPixels);

            recoveredEpsr = estimateHomogeneousEpsr(averageEpsr, ...
                backgroundEpsr, numTargetPixels, numBackgroundPixels);

            testCase.verifyEqual(recoveredEpsr, knownTargetEpsr, ...
                AbsTol=1e-12);
        end

        function testSelectionRejectsSelfNonintersectionAndInvalidDelay(testCase)
            antennaPositions = [5 21; 37 21; 5 5];
            targetMask = false(41, 41);
            targetMask(20:22, 20:22) = true;
            averageEpsr = nan(3, 3);
            averageEpsr(1, 1) = 2.0;
            averageEpsr(1, 2) = 2.0;
            averageEpsr(1, 3) = 2.0;
            pairEligible = false(3, 3);
            pairEligible(1, 1) = true;
            pairEligible(1, 2) = true;
            pairEligible(2, 1) = true;
            pairEligible(1, 3) = true;
            expectedAccepted = false(3, 3);
            expectedAccepted(1, 2) = true;

            [~, ~, ~, pairAccepted, coverage] = evaluateTargetRays( ...
                averageEpsr, 1.0, antennaPositions, targetMask, ...
                3, pairEligible);

            testCase.verifyEqual(pairAccepted, expectedAccepted);
            testCase.verifyGreaterThan(nnz(coverage), 0);
        end

        function testAngleIneligibleIntersectingPairIsRejected(testCase)
            antennaPositions = [5 21; 37 21];
            targetMask = false(41, 41);
            targetMask(20:22, 20:22) = true;
            averageEpsr = [NaN 2.0; 2.0 NaN];
            pairEligible = false(2, 2);

            [rayEstimates, ~, ~, pairAccepted, coverage] = ...
                evaluateTargetRays(averageEpsr, 1.0, antennaPositions, ...
                targetMask, 3, pairEligible);

            testCase.verifyFalse(any(pairAccepted, 'all'));
            testCase.verifyTrue(all(isnan(rayEstimates), 'all'));
            testCase.verifyEqual(nnz(coverage), 0);
        end

        function testReciprocalDirectionsRemainSeparate(testCase)
            antennaPositions = [5 21; 37 21];
            targetMask = false(41, 41);
            targetMask(20:22, 20:22) = true;
            averageEpsr = [NaN 1.5; 1.7 NaN];
            pairEligible = logical([0 1; 1 0]);

            [rayEstimates, ~, ~, pairAccepted] = evaluateTargetRays( ...
                averageEpsr, 1.0, antennaPositions, targetMask, ...
                3, pairEligible);

            testCase.verifyEqual(nnz(pairAccepted), 2);
            testCase.verifyGreaterThan(abs(rayEstimates(1, 2) - ...
                rayEstimates(2, 1)), 1e-12);
        end

        function testReconstructionUsesOnlyBackgroundAndRecoveredValue(testCase)
            targetMask = false(7, 9);
            targetMask(3:5, 4:6) = true;

            reconstruction = reconstructHomogeneousTargetMap( ...
                targetMask, 1.25, 4.75);

            testCase.verifyEqual(reconstruction(targetMask), ...
                4.75*ones(nnz(targetMask), 1), AbsTol=1e-12);
            testCase.verifyEqual(reconstruction(~targetMask), ...
                1.25*ones(nnz(~targetMask), 1), AbsTol=1e-12);
        end
    end
end

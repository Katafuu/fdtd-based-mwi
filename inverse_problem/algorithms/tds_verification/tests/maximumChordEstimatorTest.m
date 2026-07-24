classdef maximumChordEstimatorTest < matlab.unittest.TestCase
    %maximumChordEstimatorTest Tests mask chords and index-domain inversion.

    methods (TestClassSetup)
        function addVerificationLibrary(testCase)
            testDirectory = fileparts(mfilename('fullpath'));
            libraryDirectory = fullfile(fileparts(testDirectory), 'lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                libraryDirectory));
        end
    end

    methods (Test)
        function testFullMaskChordEqualsAntennaDistance(testCase)
            targetMask = true(41, 41);
            antennaPositions = [10 21; 32 21];
            pairEligible = logical([0 1; 1 0]);

            chordLengths = computeCenterlineIntersectionLengths( ...
                targetMask, antennaPositions, [1e-3 1e-3], ...
                pairEligible, 20);

            testCase.verifyEqual(chordLengths(1, 2), 22e-3, ...
                AbsTol=1e-12);
            testCase.verifyEqual(chordLengths(2, 1), 22e-3, ...
                AbsTol=1e-12);
        end

        function testSyntheticMatchedDelayRecoversPermittivity(testCase)
            measurement = maximumChordEstimatorTest.syntheticMeasurement(4.0);

            recoveredEpsr = estimateMaximumChordEpsr(measurement);

            testCase.verifyEqual(recoveredEpsr, 4.0, AbsTol=2e-3);
        end

        function testLongestChordSelectionRejectsShorterRay(testCase)
            measurement = maximumChordEstimatorTest.syntheticMeasurement(4.0);

            [~, details] = estimateMaximumChordEpsr(measurement);

            testCase.verifyEqual(details.numSelectedDirectedRays, 2);
            testCase.verifyTrue(details.selectedPairs(1, 2));
            testCase.verifyTrue(details.selectedPairs(2, 1));
            testCase.verifyFalse(details.selectedPairs(1, 3));
        end
    end

    methods (Static, Access=private)
        function measurement = syntheticMeasurement(targetEpsr)
            targetMask = false(41, 41);
            targetMask(14:28, 18:24) = true;
            antennaPositions = [5 21; 37 21; 21 5];
            pairEligible = logical([0 1 1; 1 0 0; 1 0 0]);
            chordLengths = computeCenterlineIntersectionLengths( ...
                targetMask, antennaPositions, [1e-3 1e-3], ...
                pairEligible, 40);
            dt = 1e-12;
            c0 = 3e8;
            indexContrast = sqrt(targetEpsr)-1;
            matchedDelaySteps = ...
                chordLengths .* indexContrast ./ (c0*dt);
            matchedDelaySteps(1, 3) = matchedDelaySteps(1, 3) + 100;
            matchedDelaySteps(3, 1) = matchedDelaySteps(3, 1) + 100;

            measurement = struct();
            measurement.targetMask = targetMask;
            measurement.antennaPositions = antennaPositions;
            measurement.dx = 1e-3;
            measurement.dy = 1e-3;
            measurement.pairEligible = pairEligible;
            measurement.matchedDelaySteps = matchedDelaySteps;
            measurement.c0 = c0;
            measurement.dt = dt;
            measurement.backgroundEpsr = 1.0;
        end
    end
end

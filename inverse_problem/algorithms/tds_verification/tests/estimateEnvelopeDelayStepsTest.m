classdef estimateEnvelopeDelayStepsTest < matlab.unittest.TestCase
    %estimateEnvelopeDelayStepsTest Tests normalized envelope matching.

    methods (TestClassSetup)
        function addVerificationLibrary(testCase)
            testDirectory = fileparts(mfilename('fullpath'));
            libraryDirectory = fullfile(fileparts(testDirectory), 'lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                libraryDirectory));
        end
    end

    methods (Test)
        function testIntegerDelayIsRecovered(testCase)
            time = (1:240).';
            referenceSignal = exp(-0.5*((time-80)/7).^2);
            objectSignal = exp(-0.5*((time-97)/7).^2);

            [delaySteps, score] = estimateEnvelopeDelaySteps( ...
                referenceSignal, objectSignal, 80, 24, 40);

            testCase.verifyEqual(delaySteps, 17, AbsTol=0.05);
            testCase.verifyGreaterThan(score, 0.999);
        end

        function testAmplitudeScalingDoesNotChangeDelay(testCase)
            time = (1:240).';
            referenceSignal = exp(-0.5*((time-80)/7).^2);
            objectSignal = 0.15*exp(-0.5*((time-91)/7).^2);

            delaySteps = estimateEnvelopeDelaySteps( ...
                referenceSignal, objectSignal, 80, 24, 40);

            testCase.verifyEqual(delaySteps, 11, AbsTol=0.05);
        end

        function testDelaySearchLimitIsRespected(testCase)
            time = (1:240).';
            referenceSignal = exp(-0.5*((time-80)/7).^2);
            objectSignal = exp(-0.5*((time-110)/7).^2);

            delaySteps = estimateEnvelopeDelaySteps( ...
                referenceSignal, objectSignal, 80, 24, 12);

            testCase.verifyEqual(delaySteps, 12, AbsTol=1e-12);
        end
    end
end

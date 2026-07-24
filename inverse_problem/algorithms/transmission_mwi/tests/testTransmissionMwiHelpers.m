classdef testTransmissionMwiHelpers < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addLibraryToPath(testCase)
            algorithmFolder = fileparts(fileparts(mfilename('fullpath')));
            workspaceRoot = fileparts(fileparts(fileparts(algorithmFolder)));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'buildLib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(algorithmFolder, 'lib')));
        end
    end

    methods (Test)
        function circularGeometryPreservesPairConventions(testCase)
            cfg = struct();
            cfg.grid = fdtdmat.createGrid([31 31], [1 1]);
            cfg.antennas.numAntennas = 8;
            cfg.antennas.center = [16 16];
            cfg.antennas.radius = 11;
            cfg.antennas.focusPadding = 0;
            cfg = buildCircularAntennaArrayIdx(cfg);
            antennaIdx = cfg.antennas.pos;
            center = cfg.antennas.center;

            [distance, globalAngle, betaTx, betaRx, pairWeight] = ...
                buildCircularPairGeometry(antennaIdx, 0.01, 0.01, center);

            testCase.verifySize(antennaIdx, [8 2]);
            testCase.verifyEqual(size(unique(antennaIdx, 'rows'), 1), 8);
            testCase.verifySize(distance, [8 8]);
            testCase.verifySize(globalAngle, [8 8]);
            testCase.verifySize(betaTx, [8 8]);
            testCase.verifySize(betaRx, [8 8]);
            testCase.verifySize(pairWeight, [8 8]);
            testCase.verifyTrue(all(isnan(diag(distance))));
            testCase.verifyTrue(all(isnan(diag(globalAngle))));
            testCase.verifyTrue(all(isnan(diag(betaTx))));
            testCase.verifyTrue(all(isnan(diag(betaRx))));
            testCase.verifyEqual(diag(pairWeight), zeros(8, 1));
        end

        function firstProminentPeakPrecedesLargerReflection(testCase)
            [cfg, receiverIndex, time] = ...
                testTransmissionMwiHelpers.arrivalTestConfiguration();
            precursor = testTransmissionMwiHelpers.pulse(time, 6, 1.5, 0.02);
            directPulse = testTransmissionMwiHelpers.pulse(time, 40, 4, 1);
            reflection = testTransmissionMwiHelpers.pulse(time, 100, 4, 4);
            Ez = testTransmissionMwiHelpers.receiverField(receiverIndex, ...
                precursor + directPulse + reflection, numel(time));

            arrivalStep = getArrivalSteps(Ez, receiverIndex, cfg);

            testCase.verifyEqual(arrivalStep, 40, AbsTol=1);
        end

        function signalWithoutPeakReturnsNaN(testCase)
            [cfg, receiverIndex, time] = ...
                testTransmissionMwiHelpers.arrivalTestConfiguration();
            Ez = testTransmissionMwiHelpers.receiverField( ...
                receiverIndex, zeros(size(time)), numel(time));

            arrivalStep = getArrivalSteps(Ez, receiverIndex, cfg);

            testCase.verifyTrue(isnan(arrivalStep));
        end

        function objectIncidentDifferencePreservesStepDelay(testCase)
            [cfg, receiverIndex, time] = ...
                testTransmissionMwiHelpers.arrivalTestConfiguration();
            incidentTrace = testTransmissionMwiHelpers.pulse(time, 40, 4, 1) + ...
                testTransmissionMwiHelpers.pulse(time, 100, 4, 3);
            objectTrace = testTransmissionMwiHelpers.pulse(time, 47, 4, 1) + ...
                testTransmissionMwiHelpers.pulse(time, 107, 4, 3);
            incidentEz = testTransmissionMwiHelpers.receiverField( ...
                receiverIndex, incidentTrace, numel(time));
            objectEz = testTransmissionMwiHelpers.receiverField( ...
                receiverIndex, objectTrace, numel(time));

            incidentArrival = getArrivalSteps(incidentEz, receiverIndex, cfg);
            objectArrival = getArrivalSteps(objectEz, receiverIndex, cfg);

            testCase.verifyEqual(objectArrival - incidentArrival, 7, AbsTol=1);
        end

        function multipleActiveTransmittersAreRejected(testCase)
            [cfg, receiverIndex, time] = ...
                testTransmissionMwiHelpers.arrivalTestConfiguration();
            cfg.antennas.pos = [1 1; 2 1];
            cfg.source.samples = repmat(cfg.source.samples, 2, 1);
            Ez = testTransmissionMwiHelpers.receiverField( ...
                receiverIndex, zeros(size(time)), numel(time));

            testCase.verifyError( ...
                @() getArrivalSteps(Ez, receiverIndex, cfg), ...
                'getArrivalSteps:ActiveTransmitterCount');
        end
    end

    methods (Static, Access=private)
        function [cfg, receiverIndex, time] = arrivalTestConfiguration()
            time = (1:180).';
            cfg = struct();
            cfg.c0 = 1;
            cfg.dt = 1;
            cfg.dx = 1;
            cfg.dy = 1;
            cfg.antennas.pos = [1 1];
            cfg.source.samples = ...
                testTransmissionMwiHelpers.pulse(time, 10, 3, 1).';
            receiverIndex = [11 1];
        end

        function values = pulse(time, center, width, amplitude)
            carrier = cos(0.8*(time - center));
            envelope = exp(-0.5*((time - center)/width).^2);
            values = amplitude*carrier.*envelope;
        end

        function Ez = receiverField(receiverIndex, trace, numSteps)
            Ez = zeros(receiverIndex(1), receiverIndex(2), numSteps);
            Ez(receiverIndex(1), receiverIndex(2), :) = ...
                reshape(trace, 1, 1, []);
        end
    end
end

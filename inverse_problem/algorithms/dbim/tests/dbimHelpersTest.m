classdef dbimHelpersTest < matlab.unittest.TestCase
    %dbimHelpersTest Unit tests for the FDTD-DBIM numerical helpers.

    methods (TestClassSetup)
        function addDbimPaths(testCase)
            algorithmFolder = fileparts(fileparts(mfilename('fullpath')));
            workspaceRoot = fileparts(fileparts(fileparts(algorithmFolder)));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(workspaceRoot, 'buildLib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(algorithmFolder, 'lib')));
        end
    end

    methods (Test, TestTags = {'Unit'})
        function greenFunctionIsReciprocalAndFinite(testCase)
            observation = [0 0; 0.2 -0.1];
            source = [0.1 0.3; -0.2 0.4; 0.5 0.1];
            wavenumber = 15 - 0.2i;

            forwardGreen = green2D(observation, source, wavenumber);
            reverseGreen = green2D(source, observation, wavenumber);
            coincidentGreen = green2D([0 0], [0 0], wavenumber);

            testCase.verifySize(forwardGreen, [2 3]);
            testCase.verifyEqual(forwardGreen, reverseGreen.', AbsTol=1e-12);
            testCase.verifyTrue(all(isfinite(forwardGreen(:))));
            testCase.verifyTrue(all(isfinite(coincidentGreen(:))));
        end

        function directDftRecoversComplexTone(testCase)
            sampleInterval = 0.01;
            sampleTimes = 0:sampleInterval:(1 - sampleInterval);
            frequency = 5;
            samples = 2 .* exp(1i*2*pi*frequency.*sampleTimes);

            component = extractFrequencyComponents( ...
                samples, sampleTimes, frequency);

            testCase.verifyEqual(component, 2, AbsTol=1e-12);
        end

        function directDftRejectsNonuniformTimes(testCase)
            samples = ones(1, 4);
            sampleTimes = [0 1 2 4];

            testCase.verifyError(@() extractFrequencyComponents( ...
                samples, sampleTimes, 0.1), ...
                'extractFrequencyComponents:NonuniformSampleTimes');
        end

        function jointLinearSystemMatchesComplexEquation(testCase)
            cfg = struct('mu0', 4*pi*1e-7, 'eps0', 8.854e-12, ...
                'dx', 1e-3, 'dy', 1e-3);
            frequency = 1e9;
            omega = 2*pi*frequency;
            pairs = [1 2];
            doiFields = complex(zeros(2, 2, 1));
            doiFields(:, :, 1) = [1+0.2i, 8-3i; -0.2+0.1i, -4+5i];
            analyticGreen = complex(zeros(2, 2, 1));
            analyticGreen(:, :, 1) = [0.2+0.7i, 0.6-0.1i; ...
                -0.5+0.3i, -0.4+0.8i];
            predicted = complex(zeros(2, 2, 1));
            deltaEpsr = [0.2; -0.1];
            deltaSigma = [0.01; 0.02];
            complexOperator = omega^2*cfg.mu0*cfg.eps0*cfg.dx*cfg.dy .* ...
                (analyticGreen(:, 2, 1).'.*doiFields(:, 1, 1).');
            residual = complexOperator * ...
                (deltaEpsr - 1i.*deltaSigma./(omega*cfg.eps0));
            measured = complex(zeros(2, 2, 1));
            measured(2, 1, 1) = residual;
            options = struct('operatorPrecision', "double", ...
                'epsrScale', 1, 'sigmaScale', 1, 'assemblyBlockSize', 1);
            expectedOperator = [ ...
                real(complexOperator), imag(complexOperator)./(omega*cfg.eps0); ...
                imag(complexOperator), -real(complexOperator)./(omega*cfg.eps0)];

            [operator, rhs, metadata] = buildDbimLinearSystem( ...
                measured, predicted, doiFields, pairs, analyticGreen, ...
                frequency, cfg, options);

            testCase.verifyEqual(operator, expectedOperator, AbsTol=1e-12);
            testCase.verifyEqual(rhs, [real(residual); imag(residual)], ...
                AbsTol=1e-12);
            testCase.verifyEqual(metadata.normalizedResidual, 1, AbsTol=1e-12);
            testCase.verifyEqual(metadata.greenFunction, ...
                "analytical homogeneous 2-D");
        end

        function fullTsvdRecoversNonsingularSystem(testCase)
            operator = diag([4 2 1]);
            expected = [1; -2; 0.5];
            rhs = operator * expected;
            options = struct('energyThreshold', 1, 'initialRank', 1, ...
                'maximumRank', 3, 'fullSvdMaxDimension', 10);

            [solution, info] = solve_tsvd(operator, rhs, options);

            testCase.verifyEqual(solution, expected, AbsTol=1e-12);
            testCase.verifyEqual(info.retainedRank, 3);
            testCase.verifyEqual(info.achievedEnergy, 1, AbsTol=1e-12);
            testCase.verifyFalse(info.hitRankCap);
        end

        function tsvdReportsRankCap(testCase)
            operator = diag([4 2 1]);
            rhs = operator * [1; 1; 1];
            options = struct('energyThreshold', 0.999, 'initialRank', 1, ...
                'maximumRank', 1, 'fullSvdMaxDimension', 10);

            [~, info] = solve_tsvd(operator, rhs, options);

            testCase.verifyEqual(info.retainedRank, 1);
            testCase.verifyEqual(info.achievedEnergy, 16/21, AbsTol=1e-12);
            testCase.verifyTrue(info.hitRankCap);
        end

        function constrainedUpdatePreservesExterior(testCase)
            cfg = struct();
            cfg.grid = struct( ...
                'epsr', 2.*ones(3, 3), ...
                'cond_e', 0.2.*ones(3, 3), ...
                'background', struct('epsr', 1, 'cond_e', 0));
            doiMask = false(3, 3);
            doiMask(2, 2:3) = true;
            options = struct('relaxation', 0.5, ...
                'epsrBounds', [1 3], 'sigmaBounds', [0 0.5], ...
                'epsrScale', 1, 'sigmaScale', 0.1);

            [updated, info] = applyDbimUpdate( ...
                cfg, [4; -4], [2; -2], doiMask, options);

            testCase.verifyEqual(updated.grid.epsr(doiMask), [3; 1], AbsTol=0);
            testCase.verifyEqual(updated.grid.cond_e(doiMask), [0.5; 0], AbsTol=0);
            testCase.verifyEqual(updated.grid.epsr(~doiMask), ones(7, 1), AbsTol=0);
            testCase.verifyEqual(updated.grid.cond_e(~doiMask), zeros(7, 1), AbsTol=0);
            testCase.verifyGreaterThan(info.relativeUpdateNorm, 0);
        end

    end
end

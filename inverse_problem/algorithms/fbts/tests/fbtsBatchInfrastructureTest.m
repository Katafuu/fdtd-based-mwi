classdef fbtsBatchInfrastructureTest < matlab.unittest.TestCase
    %fbtsBatchInfrastructureTest Test comparison-run option plumbing.

    properties
        TemporaryDirectory
    end

    methods (TestMethodSetup)
        function createTemporaryDirectory(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryDirectory = fixture.Folder;
        end
    end

    methods (TestClassSetup)
        function addFbtsPath(testCase)
            testFile = mfilename('fullpath');
            fbtsDirectory = fileparts(fileparts(testFile));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fbtsDirectory));
        end
    end

    methods (Test)
        function defaultAndComparisonOptionsResolve(testCase)
            defaults = resolveFbtsOptions();
            testCase.verifyEqual(defaults.numIterations, 15);
            testCase.verifyEqual( ...
                defaults.sensitivityDownsampleFactor, 4);
            testCase.verifyEqual(defaults.outputDirectory, "");

            requested = resolveFbtsOptions(struct( ...
                'numIterations', 3, ...
                'sensitivityDownsampleFactor', 1, ...
                'outputDirectory', "comparison/full", ...
                'runLabel', "without coarsening"));
            testCase.verifyEqual(requested.numIterations, 3);
            testCase.verifyEqual( ...
                requested.sensitivityDownsampleFactor, 1);
            testCase.verifyEqual( ...
                requested.outputDirectory, "comparison/full");
            testCase.verifyEqual(requested.runLabel, "without coarsening");
        end

        function unknownOptionIsRejected(testCase)
            testCase.verifyError( ...
                @() resolveFbtsOptions(struct('notAnOption', 1)), ...
                'fbts:UnknownOption');
        end

        function indexedDirectoriesUseIndependentPrefixes(testCase)
            mkdir(fullfile(testCase.TemporaryDirectory, 'run_0002'));
            mkdir(fullfile(testCase.TemporaryDirectory, ...
                'comparative_run_0004'));
            mkdir(fullfile(testCase.TemporaryDirectory, ...
                'comparative_run_notes'));

            normalDirectory = createNextFbtsOutputDirectory( ...
                testCase.TemporaryDirectory, 'run');
            comparisonDirectory = createNextFbtsOutputDirectory( ...
                testCase.TemporaryDirectory, 'comparative_run');

            testCase.verifyEqual(string(normalDirectory), string(fullfile( ...
                testCase.TemporaryDirectory, 'run_0003')));
            testCase.verifyEqual(string(comparisonDirectory), ...
                string(fullfile(testCase.TemporaryDirectory, ...
                'comparative_run_0005')));
            testCase.verifyTrue(isfolder(normalDirectory));
            testCase.verifyTrue(isfolder(comparisonDirectory));
        end
    end
end

classdef fbtsRefactorTest < matlab.unittest.TestCase
    %fbtsRefactorTest Verify function/script compatibility and saved artifacts.

    properties
        TemporaryDirectory
        FbtsDirectory
    end

    methods (TestClassSetup)
        function addRepositoryPaths(testCase)
            fbtsDir = fileparts(fileparts(mfilename('fullpath')));
            root = fileparts(fileparts(fileparts(fbtsDir)));
            testCase.FbtsDirectory = fbtsDir;
            paths = {fbtsDir, fullfile(fbtsDir, 'tests', 'helpers'), ...
                fullfile(root, 'buildLib'), ...
                fullfile(root, 'forward_solver', 'mex'), ...
                fullfile(fileparts(fbtsDir), 'time_reversal', 'lib')};
            for k = 1:numel(paths)
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(paths{k}));
            end
        end
    end

    methods (TestMethodSetup)
        function isolateOutputs(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryDirectory = fixture.Folder;
            previousVisibility = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            testCase.addTeardown(@() set(groot, ...
                'DefaultFigureVisible', previousVisibility));
        end
    end

    methods (Test)
        function optionsKeepDefaultsAndRejectInvalidControls(testCase)
            defaults = resolveFbtsOptions();
            testCase.verifyEqual(defaults.numIterations, 15);
            testCase.verifyEqual(defaults.sensitivityDownsampleFactor, 4);
            testCase.verifyEqual(defaults.outputDirectory, "");
            testCase.verifyEqual(defaults.runLabel, "");
            testCase.verifyEqual(resolveFbtsOptions([]), defaults);
            testCase.verifyError(@() resolveFbtsOptions(1), 'fbts:InvalidOptions');
            testCase.verifyError(@() resolveFbtsOptions(struct('unknown', 1)), ...
                'fbts:UnknownOption');
            testCase.verifyError(@() resolveFbtsOptions(struct('runLabel', {1})), ...
                'fbts:InvalidTextOption');
        end

        function functionsAndScriptSaveIdenticalResults(testCase)
            originalCfg = makeFbtsTestConfig();
            options = struct('numIterations', 3, ...
                'sensitivityDownsampleFactor', 4, 'runLabel', "regression", ...
                'outputDirectory', string(fullfile(testCase.TemporaryDirectory, 'direct')));
            beforeFigures = findall(groot, 'Type', 'figure');
            [computed, computedCfg] = runFbts(originalCfg, options);
            testCase.verifyFalse(isfolder(options.outputDirectory));
            testCase.verifyFalse(isfield(computed, 'output_files'));
            testCase.verifyEqual(findall(groot, 'Type', 'figure'), beforeFigures);
            testCase.verifyFalse(isfield(originalCfg.source, 'func'));
            testCase.verifyTrue(isfield(computedCfg.source, 'func'));

            % Saving must take the iteration count from results, not options.
            savingOptions = options;
            savingOptions.numIterations = 99;
            [saved, handles] = saveFbtsResults(computedCfg, computed, savingOptions);
            testCase.addTeardown(@() closeValidFigures(struct2array(handles)));
            testCase.verifyEqual(fieldnames(handles), ...
                {'trueFigure'; 'estimatedFigure'; 'convergenceFigure'; 'costFigure'});
            testCase.verifyTrue(all(isgraphics(struct2array(handles))));
            ax = findall(handles.estimatedFigure, 'Type', 'axes');
            % Preserve plot_fbts's existing literal title and numeric subtitle.
            testCase.verifyEqual(ax.Title.String, ...
                'FBTS relative permittivity reconstruction, %d iterations');
            testCase.verifyEqual(ax.Subtitle.String, '3');

            cfg = originalCfg;
            fbtsOptions = options;
            fbtsOptions.outputDirectory = string(fullfile(testCase.TemporaryDirectory, 'script'));
            mkdir(fbtsOptions.outputDirectory);
            fbts_demo;
            testCase.addTeardown(@() closeValidFigures( ...
                [trueFigure estimatedFigure convergenceFigure costFigure]));
            testCase.verifyEqual(activeFbtsOptions, resolveFbtsOptions(fbtsOptions));
            testCase.verifyEqual(numIterations, 3);
            ignored = {'iteration_runtime', 'total_runtime', 'average_runtime', ...
                'output_directory', 'output_files', 'timing'};
            testCase.verifyEqual(rmfield(results, ignored), rmfield(saved, ignored));
            testCase.verifyEqual(func2str(cfg.source.func), func2str(computedCfg.source.func));
            testCase.verifyEqual(cfg.source.func((0:cfg.Nt-1)*cfg.dt), ...
                computedCfg.source.func((0:cfg.Nt-1)*cfg.dt));
            testCase.verifyEqual(rmfield(cfg, 'source'), rmfield(computedCfg, 'source'));
            testCase.verifyEqual(results.total_runtime, sum(results.iteration_runtime));
            testCase.verifyEqual(results.average_runtime, mean(results.iteration_runtime));
            testCase.verifyTrue(all(isfinite(results.iteration_runtime)));
            testCase.verifyGreaterThanOrEqual(results.iteration_runtime, zeros(3, 1));

            stored = load(results.output_files.mat_file);
            testCase.verifyEqual(sort(fieldnames(stored)), {'cfg'; 'results'});
            testCase.verifyEqual(stored.results, results);
            testCase.verifyEqual(rmfield(stored.cfg, 'source'), rmfield(cfg, 'source'));
            testCase.verifyEqual(stored.cfg.source.func((0:cfg.Nt-1)*cfg.dt), ...
                cfg.source.func((0:cfg.Nt-1)*cfg.dt));
            for k = 1:4
                testCase.verifyEqual(imread(results.output_files.figures(k)), ...
                    imread(saved.output_files.figures(k)));
            end
        end

        function defaultOutputUsesNextDirectory(testCase)
            % Isolate the default output root without touching existing runs.
            copyfile(fullfile(testCase.FbtsDirectory, 'saveFbtsResults.m'), ...
                testCase.TemporaryDirectory);
            copyfile(fullfile(testCase.FbtsDirectory, 'fbts_demo.m'), ...
                testCase.TemporaryDirectory);
            previousDirectory = pwd;
            testCase.addTeardown(@() cd(previousDirectory));
            cd(testCase.TemporaryDirectory);
            mkdir(fullfile(testCase.TemporaryDirectory, 'figs', 'run_0002'));
            cfg = makeFbtsTestConfig();
            [computed, ~] = runFbts(cfg);
            fbts_demo;
            testCase.addTeardown(@() closeValidFigures( ...
                [trueFigure estimatedFigure convergenceFigure costFigure]));
            testCase.verifyEqual(numIterations, 15);
            testCase.verifyEqual(activeFbtsOptions, resolveFbtsOptions());
            testCase.verifyEqual(results.epsr_est, computed.epsr_est);
            testCase.verifyEqual(results.output_directory, string(fullfile( ...
                testCase.TemporaryDirectory, 'figs', 'run_0003')));
            [saved, handles] = saveFbtsResults(cfg, results);
            testCase.addTeardown(@() closeValidFigures(struct2array(handles)));
            testCase.verifyEqual(saved.output_directory, string(fullfile( ...
                testCase.TemporaryDirectory, 'figs', 'run_0004')));
            testCase.verifyTrue(isfile(saved.output_files.mat_file));
        end

        function nonemptyDirectoryIsRejectedBeforePlotting(testCase)
            mkdir(fullfile(testCase.TemporaryDirectory, 'existing'));
            testCase.verifyError(@() saveFbtsResults(struct(), ...
                struct('num_iterations', 1), ...
                struct('outputDirectory', testCase.TemporaryDirectory)), ...
                'fbts:OutputDirectoryNotEmpty');
        end

        function missingConfigurationKeepsOriginalError(testCase)
            testCase.verifyError(@() runFbts(), 'fbts:MissingConfig');
            testCase.verifyError(@runScriptWithoutConfig, 'fbts:MissingConfig');
        end
    end
end

function runScriptWithoutConfig()
fbts_demo;
end

function closeValidFigures(handles)
close(handles(isgraphics(handles)));
end

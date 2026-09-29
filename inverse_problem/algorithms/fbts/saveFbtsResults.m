function [results, figureHandles] = saveFbtsResults(cfg, results, fbtsOptions)
%saveFbtsResults Plot and save cfg/results using the existing FBTS layout.
%   Optional fbtsOptions.outputDirectory selects an empty output directory.
%   Otherwise the next figs/run_XXXX directory is created beside this file.
%   Figures remain open and their named handles are returned to the caller.

if nargin < 3
    fbtsOptions = [];
end
activeFbtsOptions = resolveFbtsOptions(fbtsOptions);
numIterations = results.num_iterations; %#ok<NASGU> Used by plot_fbts.

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end
runOutputDirectory = createFbtsOutputDirectory( ...
    activeFbtsOptions.outputDirectory, fullfile(scriptDir, 'figs'), 'run');
figureFiles = [ ...
    string(fullfile(runOutputDirectory, ...
        '01_true_relative_permittivity.png'))
    string(fullfile(runOutputDirectory, ...
        '02_estimated_relative_permittivity.png'))
    string(fullfile(runOutputDirectory, ...
        '03_relative_error_convergence.png'))
    string(fullfile(runOutputDirectory, ...
        '04_cost_convergence.png'))
];
resultFile = string(fullfile(runOutputDirectory, 'fbts_run_data.mat'));
results.output_directory = string(runOutputDirectory);
results.output_files = struct( ...
    'figures', figureFiles, ...
    'mat_file', resultFile);

save(char(resultFile), 'cfg', 'results', '-v7.3');
plot_fbts;
fprintf('Saved FBTS figures and run data under %s.\n', ...
    runOutputDirectory);

figureHandles = struct( ...
    'trueFigure', trueFigure, ...
    'estimatedFigure', estimatedFigure, ...
    'convergenceFigure', convergenceFigure, ...
    'costFigure', costFigure);
end

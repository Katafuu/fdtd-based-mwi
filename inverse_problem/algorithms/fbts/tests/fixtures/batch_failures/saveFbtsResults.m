function [results, handles] = saveFbtsResults(cfg, results, options)
% Test double: leave an open figure and fail before returning its handle.
mkdir(options.outputDirectory);
handle = figure('Visible', 'off');
if isequal(cfg.antennas.originalIndices, 2)
    error('fbtsTest:SaveFailure', 'Injected plotting failure.');
end
if isequal(cfg.antennas.originalIndices, [1 3])
    mkdir(fullfile(options.outputDirectory, 'failure.mat'));
    error('fbtsTest:DiagnosticFailure', 'Injected diagnostic-file conflict.');
end
save(fullfile(options.outputDirectory, 'fbts_run_data.mat'), 'cfg', 'results');
handles = struct('trueFigure', handle);
end

% main Run one FBTS reconstruction and save its MATLAB data.
mainTimer = tic;
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'build_cfg.m'));
[results, cfg] = runFbts(cfg);

runOutputDirectory = createNextRunDirectory(fullfile(scriptDir, 'figs'));
resultFile = string(fullfile(runOutputDirectory, 'fbts_run_data.mat'));
results.output_directory = string(runOutputDirectory);
results.output_files = struct('mat_file', resultFile);

writeTimer = tic;
save_results(resultFile, cfg, results);
writeTime = toc(writeTimer);
fileInfo = dir(char(resultFile));
totalExecutionTime = toc(mainTimer);
fprintf('Saved FBTS run data under %s.\n', runOutputDirectory);
fprintf('Total execution time: %.3f seconds.\n', totalExecutionTime);
fprintf('MAT write time: %.3f seconds.\n', writeTime);
fprintf('MAT file size: %d bytes (%.3f MiB).\n', ...
    fileInfo.bytes, fileInfo.bytes / 2^20);

function runDirectory = createNextRunDirectory(outputRoot)
outputRoot = char(outputRoot);
if ~isfolder(outputRoot)
    [created, message] = mkdir(outputRoot);
    if ~created
        error('fbts:CreateOutputRootFailed', '%s', message);
    end
end

existingRuns = dir(fullfile(outputRoot, 'run_*'));
existingRuns = existingRuns([existingRuns.isdir]);
runNumbers = nan(numel(existingRuns), 1);
numMatchedRuns = 0;
for directoryIndex = 1:numel(existingRuns)
    token = regexp(existingRuns(directoryIndex).name, ...
        '^run_(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        numMatchedRuns = numMatchedRuns + 1;
        runNumbers(numMatchedRuns) = str2double(token{1});
    end
end
runNumbers = runNumbers(1:numMatchedRuns);
if isempty(runNumbers)
    nextRunNumber = 1;
else
    nextRunNumber = max(runNumbers) + 1;
end

runDirectory = fullfile(outputRoot, sprintf('run_%04d', nextRunNumber));
[created, message] = mkdir(runDirectory);
if ~created
    error('fbts:CreateRunDirectoryFailed', '%s', message);
end
end

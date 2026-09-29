function outputDirectory = createFbtsOutputDirectory(requestedDirectory, outputRoot, prefix)
%createFbtsOutputDirectory Use an empty explicit folder or the next prefix_XXXX.

if strlength(requestedDirectory) > 0
    outputDirectory = char(requestedDirectory);
    if isfolder(outputDirectory)
        existingOutput = dir(outputDirectory);
        existingOutput = existingOutput(~ismember( ...
            {existingOutput.name}, {'.', '..'}));
        if ~isempty(existingOutput)
            error('fbts:OutputDirectoryNotEmpty', ...
                'The requested output directory is not empty: %s', outputDirectory);
        end
    else
        [created, message] = mkdir(outputDirectory);
        if ~created
            error('fbts:CreateRunDirectoryFailed', '%s', message);
        end
    end
    return
end

if ~isfolder(outputRoot)
    [created, message] = mkdir(outputRoot);
    if ~created
        error('fbts:CreateOutputRootFailed', '%s', message);
    end
end
existingRuns = dir(fullfile(outputRoot, [prefix '_*']));
existingRuns = existingRuns([existingRuns.isdir]);
nextRunNumber = 1;
expression = ['^' regexptranslate('escape', prefix) '_(\d+)$'];
for directoryIndex = 1:numel(existingRuns)
    token = regexp(existingRuns(directoryIndex).name, expression, 'tokens', 'once');
    if ~isempty(token)
        nextRunNumber = max(nextRunNumber, str2double(token{1}) + 1);
    end
end
outputDirectory = fullfile(outputRoot, sprintf('%s_%04d', prefix, nextRunNumber));
[created, message] = mkdir(outputDirectory);
if ~created
    error('fbts:CreateRunDirectoryFailed', '%s', message);
end
end

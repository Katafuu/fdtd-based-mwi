function runDirectory = createNextFbtsOutputDirectory(outputRoot, prefix)
%createNextFbtsOutputDirectory Create the next indexed FBTS output folder.

outputRoot = char(string(outputRoot));
prefix = char(string(prefix));
if isempty(outputRoot)
    error('fbts:EmptyOutputRoot', 'The output root cannot be empty.');
end
if isempty(prefix) || contains(prefix, '/') || contains(prefix, '\')
    error('fbts:InvalidOutputPrefix', ...
        'The output directory prefix must be a nonempty folder name.');
end
if ~isfolder(outputRoot)
    [created, message] = mkdir(outputRoot);
    if ~created
        error('fbts:CreateOutputRootFailed', '%s', message);
    end
end

existingRuns = dir(fullfile(outputRoot, [prefix '_*']));
existingRuns = existingRuns([existingRuns.isdir]);
runNumbers = nan(numel(existingRuns), 1);
numMatchedRuns = 0;
expression = ['^' regexptranslate('escape', prefix) '_(\d+)$'];
for directoryIndex = 1:numel(existingRuns)
    token = regexp(existingRuns(directoryIndex).name, ...
        expression, 'tokens', 'once');
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

runDirectory = fullfile(outputRoot, ...
    sprintf('%s_%04d', prefix, nextRunNumber));
[created, message] = mkdir(runDirectory);
if ~created
    error('fbts:CreateRunDirectoryFailed', '%s', message);
end
end

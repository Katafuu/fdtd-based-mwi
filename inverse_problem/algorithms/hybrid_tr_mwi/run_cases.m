% run_cases Build one shared cfg and run both headless focus-bias cases.

controllerDir = fileparts(mfilename('fullpath'));
addpath(controllerDir, '-begin');
cfg = build_cfg(struct());

caseRoots = {fullfile(controllerDir, 'figs')};
existingRunNumbers = zeros(0, 1);

for rootIndex = 1:numel(caseRoots)
    if ~isfolder(caseRoots{rootIndex})
        continue;
    end

    runDirectories = dir(fullfile(caseRoots{rootIndex}, 'run_*'));
    runDirectories = runDirectories([runDirectories.isdir]);
    for directoryIndex = 1:numel(runDirectories)
        token = regexp(runDirectories(directoryIndex).name, ...
            '^run_(\d+)$', 'tokens', 'once');
        if ~isempty(token)
            existingRunNumbers(end + 1, 1) = str2double(token{1}); %#ok<SAGROW>
        end
    end
end

if isempty(existingRunNumbers)
    nextRunNumber = 1;
else
    nextRunNumber = max(existingRunNumbers) + 1;
end
caseRunName = sprintf('run_%04d', nextRunNumber);

fprintf('Running original and TR-enhanced cases as %s.\n', caseRunName);
run(fullfile(controllerDir, 'tds_2d_original.m'));
run(fullfile(controllerDir, 'tds_2d_tr_enhanced.m'));
fprintf('Saved both cases under figs using %s.\n', caseRunName);

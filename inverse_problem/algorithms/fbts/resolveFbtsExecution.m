function execution = resolveFbtsExecution(candidate)
%resolveFbtsExecution Optional iteration I/O; separate from numerical options.
execution = struct('checkpointFile', '', 'checkpointEvery', 1, ...
    'historyFile', '', 'resume', false, 'onIteration', [], 'computeImageError', true, ...
    'expectedCodeHash', '');
if nargin == 0 || isempty(candidate), return; end
assert(isstruct(candidate) && isscalar(candidate), 'fbts:InvalidExecution', 'Expected scalar execution struct.');
names = fieldnames(candidate);
for k = 1:numel(names)
    assert(isfield(execution,names{k}), 'fbts:InvalidExecution', 'Unknown execution option %s.',names{k});
    execution.(names{k}) = candidate.(names{k});
end
validateattributes(execution.checkpointEvery, {'numeric'}, {'scalar','integer','positive'});
for name = {'checkpointFile','historyFile','expectedCodeHash'}
    value = execution.(name{1});
    assert((ischar(value) && (isrow(value) || isempty(value))) || ...
        (isstring(value) && isscalar(value) && ~ismissing(value)), ...
        'fbts:InvalidExecution','Paths must be text scalars.');
    execution.(name{1}) = char(value);
end
for name = {'resume','computeImageError'}
    validateattributes(execution.(name{1}), {'logical'}, {'scalar'});
end
assert(isempty(execution.onIteration) || isa(execution.onIteration,'function_handle'), ...
    'fbts:InvalidExecution','onIteration must be a function handle.');
if ~isempty(execution.checkpointFile) && isempty(execution.historyFile)
    execution.historyFile = fullfile(fileparts(execution.checkpointFile),'history_work.mat');
end
end

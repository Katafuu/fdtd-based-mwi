function options = resolveFbtsBatchOptions(candidate)
%resolveFbtsBatchOptions Scheduling/storage controls, separate from FBTS math.
options = struct('symmetry','dihedral','seed',0,'disabledCounts',[], ...
    'maxCandidates',2e6,'resume',false,'dryRun',false,'checkpointEvery',1, ...
    'saveFigures',false,'maxWorkers',1,'maxCases',Inf);
if nargin == 0 || isempty(candidate), return; end
assert(isstruct(candidate) && isscalar(candidate),'fbts:InvalidBatchOptions','Expected a scalar struct.');
for key = fieldnames(candidate).'
    assert(isfield(options,key{1}),'fbts:InvalidBatchOptions','Unknown batch option %s.',key{1});
    options.(key{1}) = candidate.(key{1});
end
options.symmetry = validatestring(options.symmetry,{'dihedral','rotation','none'});
validateattributes(options.seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(options.maxCandidates,{'numeric'},{'scalar','positive'});
validateattributes(options.maxCases,{'numeric'},{'scalar','nonnegative'});
assert(isinf(options.maxCases) || options.maxCases==fix(options.maxCases), ...
    'fbts:InvalidBatchOptions','maxCases must be an integer or Inf.');
for key = {'maxWorkers','checkpointEvery'}
    validateattributes(options.(key{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
for key = {'resume','dryRun','saveFigures'}
    validateattributes(options.(key{1}),{'logical'},{'scalar'});
end
end

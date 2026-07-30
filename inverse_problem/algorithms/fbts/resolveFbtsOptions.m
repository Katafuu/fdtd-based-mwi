function options = resolveFbtsOptions(candidate)
%resolveFbtsOptions Validate optional controls for the FBTS script.

options = struct( ...
    'numIterations', 15, ...
    'sensitivityDownsampleFactor', 4, ...
    'outputDirectory', "", ...
    'runLabel', "");

if nargin == 0 || isempty(candidate)
    return
end
if ~isstruct(candidate) || ~isscalar(candidate)
    error('fbts:InvalidOptions', ...
        'fbtsOptions must be a scalar structure.');
end

candidateFields = fieldnames(candidate);
unknownFields = setdiff(candidateFields, fieldnames(options));
if ~isempty(unknownFields)
    error('fbts:UnknownOption', ...
        'Unknown FBTS option: %s.', unknownFields{1});
end
for fieldIndex = 1:numel(candidateFields)
    fieldName = candidateFields{fieldIndex};
    options.(fieldName) = candidate.(fieldName);
end

validateattributes(options.numIterations, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'numIterations');
validateattributes(options.sensitivityDownsampleFactor, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'sensitivityDownsampleFactor');

options.outputDirectory = validateTextScalar( ...
    options.outputDirectory, 'outputDirectory');
options.runLabel = validateTextScalar(options.runLabel, 'runLabel');
end

function value = validateTextScalar(value, argumentName)
if ~(ischar(value) && (isrow(value) || isempty(value))) && ...
        ~(isstring(value) && isscalar(value) && ~ismissing(value))
    error('fbts:InvalidTextOption', ...
        '%s must be a character vector or string scalar.', argumentName);
end
value = string(value);
end

function [cfgUpdated, info] = applyDbimUpdate( ...
    cfgCurrent, deltaEpsr, deltaSigma, doiMask, options)
%applyDbimUpdate Apply a constrained material update inside the DOI.

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(cfgCurrent) || ~isfield(cfgCurrent, 'grid') || ...
        ~all(isfield(cfgCurrent.grid, {'epsr', 'cond_e', 'background'}))
    error('applyDbimUpdate:InvalidConfig', ...
        'cfgCurrent.grid must contain epsr, cond_e, and background.');
end
doiMask = logical(doiMask);
if ~isequal(size(doiMask), size(cfgCurrent.grid.epsr)) || ...
        ~isequal(size(cfgCurrent.grid.cond_e), size(cfgCurrent.grid.epsr))
    error('applyDbimUpdate:SizeMismatch', ...
        'The DOI and material arrays must have matching sizes.');
end
numDoiPixels = nnz(doiMask);
deltaEpsr = double(deltaEpsr(:));
deltaSigma = double(deltaSigma(:));
if numel(deltaEpsr) ~= numDoiPixels || numel(deltaSigma) ~= numDoiPixels || ...
        any(~isfinite(deltaEpsr)) || any(~isfinite(deltaSigma))
    error('applyDbimUpdate:InvalidUpdate', ...
        'Update vectors must contain one finite value per DOI pixel.');
end

relaxation = scalarOption(options, 'relaxation', 0.5, 0, 1);
epsrBounds = boundsOption(options, 'epsrBounds', [1 80]);
sigmaBounds = boundsOption(options, 'sigmaBounds', [0 2]);
epsrScale = scalarOption(options, 'epsrScale', 1.0, 0, Inf);
sigmaScale = scalarOption(options, 'sigmaScale', 0.1, 0, Inf);

epsrBefore = double(cfgCurrent.grid.epsr);
sigmaBefore = double(cfgCurrent.grid.cond_e);
epsrAfter = epsrBefore;
sigmaAfter = sigmaBefore;

epsrAfter(doiMask) = min(max( ...
    epsrBefore(doiMask) + relaxation .* deltaEpsr, epsrBounds(1)), ...
    epsrBounds(2));
sigmaAfter(doiMask) = min(max( ...
    sigmaBefore(doiMask) + relaxation .* deltaSigma, sigmaBounds(1)), ...
    sigmaBounds(2));

background = cfgCurrent.grid.background;
epsrAfter(~doiMask) = background.epsr;
sigmaAfter(~doiMask) = background.cond_e;

cfgUpdated = cfgCurrent;
cfgUpdated.grid.epsr = epsrAfter;
cfgUpdated.grid.cond_e = sigmaAfter;

appliedEpsr = epsrAfter(doiMask) - epsrBefore(doiMask);
appliedSigma = sigmaAfter(doiMask) - sigmaBefore(doiMask);
scaledUpdate = [appliedEpsr ./ epsrScale; appliedSigma ./ sigmaScale];
scaledState = [epsrBefore(doiMask) ./ epsrScale; ...
    sigmaBefore(doiMask) ./ sigmaScale];
info = struct();
info.appliedDeltaEpsr = appliedEpsr;
info.appliedDeltaSigma = appliedSigma;
info.relativeUpdateNorm = norm(scaledUpdate) / max(norm(scaledState), eps);
info.numClampedEpsr = nnz(epsrAfter(doiMask) == epsrBounds(1) | ...
    epsrAfter(doiMask) == epsrBounds(2));
info.numClampedSigma = nnz(sigmaAfter(doiMask) == sigmaBounds(1) | ...
    sigmaAfter(doiMask) == sigmaBounds(2));
end

function value = scalarOption(options, name, defaultValue, lowerBound, upperBound)
value = getOption(options, name, defaultValue);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value <= lowerBound || value > upperBound
    error('applyDbimUpdate:InvalidOption', ...
        'options.%s must lie in (%g, %g].', name, lowerBound, upperBound);
end
value = double(value);
end

function value = boundsOption(options, name, defaultValue)
value = double(getOption(options, name, defaultValue));
value = value(:).';
if numel(value) ~= 2 || any(~isfinite(value)) || value(2) < value(1)
    error('applyDbimUpdate:InvalidOption', ...
        'options.%s must be a finite ordered two-element vector.', name);
end
end

function value = getOption(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

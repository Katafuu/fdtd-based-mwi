function [cfgVaried, info] = applyEpsrInhomogeneity(cfg, mask, opts)
%applyEpsrInhomogeneity Add bounded spatially correlated epsilon-r offsets.
% The offset is centered on the existing permittivity in each selected cell.
% opts.variance and opts.correlationCells are required; opts.seed is optional.
% info reports offset statistics and the resulting selected-cell epsilon-r range.
if ~isstruct(cfg) || ~isscalar(cfg) || ~isfield(cfg, 'grid') || ...
        ~isfield(cfg.grid, 'epsr') || ~isnumeric(cfg.grid.epsr) || ...
        ~isreal(cfg.grid.epsr) || any(~isfinite(cfg.grid.epsr(:)))
    error('fdtdmat:applyEpsrInhomogeneity:InvalidConfig', ...
        'cfg.grid.epsr must be a finite real numeric map.');
end
if ~(islogical(mask) || isnumeric(mask)) || ...
        ~isequal(size(mask), size(cfg.grid.epsr)) || ...
        any(~isfinite(double(mask(:))))
    error('fdtdmat:applyEpsrInhomogeneity:InvalidMask', ...
        'mask must be a finite array matching cfg.grid.epsr.');
end
mask = logical(mask);
if ~any(mask(:))
    error('fdtdmat:applyEpsrInhomogeneity:EmptyMask', ...
        'mask must select at least one cell.');
end
if ~isstruct(opts) || ~isscalar(opts) || ...
        ~isfield(opts, 'variance') || ~isfield(opts, 'correlationCells')
    error('fdtdmat:applyEpsrInhomogeneity:InvalidOptions', ...
        'opts must contain variance and correlationCells.');
end
validateattributes(opts.variance, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'nonnegative'});
validateattributes(opts.correlationCells, {'numeric'}, ...
    {'real', 'scalar', 'finite', 'positive'});
seed = [];
if isfield(opts, 'seed') && ~isempty(opts.seed)
    seed = opts.seed;
    validateattributes(seed, {'numeric'}, ...
        {'real', 'scalar', 'integer', '>=', 0, '<=', 2^32-1});
end

cfgVaried = cfg;
offset = zeros(nnz(mask), 1);
if opts.variance > 0
    if ~isempty(seed)
        priorRng = rng;
        rngCleanup = onCleanup(@() rng(priorRng)); %#ok<NASGU>
        rng(seed, 'twister');
    end
    kernelRadius = ceil(3 * opts.correlationCells);
    kernelAxis = -kernelRadius:kernelRadius;
    kernel = exp(-0.5 * (kernelAxis ./ opts.correlationCells).^2);
    kernel = kernel ./ sum(kernel);
    smoothField = conv2(kernel(:), kernel(:).', ...
        randn(size(mask)), 'same');
    maskedValues = smoothField(mask);
    maskedStd = std(maskedValues);
    if maskedStd == 0
        error('fdtdmat:applyEpsrInhomogeneity:ConstantRandomField', ...
            'The random field has no variation within the selected region.');
    end
    normalizedField = (maskedValues - mean(maskedValues)) ./ maskedStd;
    normalizedField = max(-2.5, min(2.5, normalizedField));
    priorValues = cfg.grid.epsr(mask);
    updatedValues = max(1.0, priorValues + ...
        sqrt(opts.variance) .* normalizedField);
    cfgVaried.grid.epsr(mask) = updatedValues;
    offset = updatedValues - priorValues;
end
selectedValues = cfgVaried.grid.epsr(mask);
info = struct( ...
    'seed', seed, ...
    'requestedVariance', opts.variance, ...
    'realizedVarianceOfOffset', var(offset), ...
    'meanOffset', mean(offset), ...
    'correlationCells', opts.correlationCells, ...
    'modifiedCellCount', nnz(offset ~= 0), ...
    'minEpsr', min(selectedValues), ...
    'maxEpsr', max(selectedValues));
end

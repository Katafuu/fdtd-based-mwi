function grid = applyRegion(grid, region, properties)
%applyRegion Apply material properties to cells selected by region.mask.

if ~isstruct(properties)
    error('fdtdmat:applyRegion:InvalidProperties', 'properties must be a struct.');
end

supportedProperties = {'epsr', 'murx', 'mury', 'cond_e', 'cond_m'};
propertyNames = fieldnames(properties);
for propertyIndex = 1:numel(propertyNames)
    propertyName = propertyNames{propertyIndex};
    if ~ismember(propertyName, supportedProperties)
        error('fdtdmat:applyRegion:UnsupportedProperty', ...
            'Unsupported material property "%s".', propertyName);
    end

    value = properties.(propertyName);
    grid.(propertyName) = assignMaskedValue(grid.(propertyName), region.mask, value, propertyName);
end
end

function matrix = assignMaskedValue(matrix, mask, value, propertyName)
targetMask = projectedMask(mask, size(matrix), propertyName);

if isscalar(value)
    matrix(targetMask) = value;
elseif isequal(size(value), size(matrix))
    matrix(targetMask) = value(targetMask);
else
    error('fdtdmat:applyRegion:ValueSizeMismatch', ...
        'Value for "%s" must be scalar or the same size as the material map.', propertyName);
end
end

function targetMask = projectedMask(mask, matrixSize, propertyName)
maskSize = size(mask);
if isequal(matrixSize, maskSize)
    targetMask = mask;
    return;
end

if strcmp(propertyName, 'murx') && matrixSize(1) == maskSize(1) && matrixSize(2) == maskSize(2) - 1
    targetMask = mask(:, 1:end-1);
    return;
end

if strcmp(propertyName, 'mury') && matrixSize(1) == maskSize(1) - 1 && matrixSize(2) == maskSize(2)
    targetMask = mask(1:end-1, :);
    return;
end

error('fdtdmat:applyRegion:MaskSizeMismatch', ...
    'Region mask cannot be projected onto "%s" map size.', propertyName);
end

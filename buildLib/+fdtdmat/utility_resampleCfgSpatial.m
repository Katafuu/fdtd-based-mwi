function resizedCfg = utility_resampleCfgSpatial(cfg, factor, direction)
%utility_resampleCfgSpatial Shared implementation for cfg spatial resampling.

validateInputs(cfg, factor, direction);
factor = double(factor);
direction = char(direction);
if factor == 1
    resizedCfg = cfg;
    return
end

oldNx = positiveInteger(cfg.Nx, 'fdtdmat:resampleCfg:InvalidNx');
oldNy = positiveInteger(cfg.Ny, 'fdtdmat:resampleCfg:InvalidNy');
oldDx = positiveScalar(cfg.dx, 'fdtdmat:resampleCfg:InvalidDx');
oldDy = positiveScalar(cfg.dy, 'fdtdmat:resampleCfg:InvalidDy');

if strcmp(direction, 'down')
    if mod(oldNx, factor) ~= 0 || mod(oldNy, factor) ~= 0
        error('fdtdmat:resampleCfg:NondivisibleGrid', ...
            'cfg.Nx and cfg.Ny must be divisible by the downsampling factor.');
    end
    newNx = oldNx / factor;
    newNy = oldNy / factor;
    newDx = oldDx * factor;
    newDy = oldDy * factor;
else
    newNx = oldNx * factor;
    newNy = oldNy * factor;
    newDx = oldDx / factor;
    newDy = oldDy / factor;
end
if newNx < 2 || newNy < 2
    error('fdtdmat:resampleCfg:GridTooSmall', ...
        'The resized grid must contain at least two samples on each axis.');
end

oldSpec = struct('Nx', oldNx, 'Ny', oldNy, 'dx', oldDx, 'dy', oldDy);
newSpec = struct('Nx', newNx, 'Ny', newNy, 'dx', newDx, 'dy', newDy);

resizedCfg = cfg;
resizedCfg.Nx = newNx;
resizedCfg.Ny = newNy;
resizedCfg.dx = newDx;
resizedCfg.dy = newDy;
resizedCfg.grid = resizeGrid(cfg.grid, oldSpec, newSpec, factor, direction);

if isfield(cfg, 'pml') && ~isempty(cfg.pml)
    resizedCfg.pml = resizePml(cfg.pml, oldSpec, newSpec, factor, direction);
end
if isfield(cfg, 'antennas') && ~isempty(cfg.antennas)
    resizedCfg.antennas = resizeAntennas( ...
        cfg.antennas, oldSpec, newSpec, factor, direction);
end
if isfield(cfg, 'targets') && ~isempty(cfg.targets)
    resizedCfg.targets = resizeTargets( ...
        cfg.targets, oldSpec, newSpec, factor, direction);
end
if isfield(cfg, 'source') && isstruct(cfg.source)
    resizedCfg.source = cfg.source;
    if isfield(cfg.source, 'location') && ~isempty(cfg.source.location)
        resizedCfg.source.location = mapIndexCoordinates( ...
            cfg.source.location, oldSpec, newSpec, true);
    end
    if isfield(resizedCfg, 'antennas') && ...
            isfield(resizedCfg.antennas, 'pos') && ...
            isfield(resizedCfg.source, 'location')
        resizedCfg.source.location = resizedCfg.antennas.pos;
    end
end
end

function validateInputs(cfg, factor, direction)
if ~isstruct(cfg) || ~isscalar(cfg)
    error('fdtdmat:resampleCfg:InvalidCfg', 'cfg must be a scalar struct.');
end
requiredFields = {'Nx', 'Ny', 'dx', 'dy', 'grid'};
for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields{fieldIndex};
    if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
        error('fdtdmat:resampleCfg:MissingCfgField', ...
            'cfg.%s is required.', fieldName);
    end
end
if ~isnumeric(factor) || ~isscalar(factor) || ~isfinite(factor) || ...
        factor < 1 || factor ~= round(factor)
    error('fdtdmat:resampleCfg:InvalidFactor', ...
        'factor must be a positive integer scalar.');
end
direction = char(string(direction));
if ~ismember(direction, {'down', 'up'})
    error('fdtdmat:resampleCfg:InvalidDirection', ...
        'direction must be "down" or "up".');
end
end

function newGrid = resizeGrid(oldGrid, oldSpec, newSpec, factor, direction)
if ~isstruct(oldGrid) || ~isscalar(oldGrid) || ...
        ~isfield(oldGrid, 'background')
    error('fdtdmat:resampleCfg:InvalidGrid', ...
        'cfg.grid must be a scalar grid struct with a background field.');
end
originPhysical = [0 0];
if isfield(oldGrid, 'originPhysical') && ~isempty(oldGrid.originPhysical)
    originPhysical = oldGrid.originPhysical;
end
newGrid = fdtdmat.createGrid( ...
    [newSpec.Nx newSpec.Ny], [newSpec.dx newSpec.dy], ...
    oldGrid.background, originPhysical);

cellFields = {'epsr', 'cond_e', 'cond_m', ...
    'epsr_bg', 'cond_e_bg', 'cond_m_bg'};
murxFields = {'murx', 'murx_bg'};
muryFields = {'mury', 'mury_bg'};
for fieldIndex = 1:numel(cellFields)
    newGrid = copyGridMap(newGrid, oldGrid, cellFields{fieldIndex}, ...
        'cell', oldSpec, newSpec, factor, direction);
end
for fieldIndex = 1:numel(murxFields)
    newGrid = copyGridMap(newGrid, oldGrid, murxFields{fieldIndex}, ...
        'murx', oldSpec, newSpec, factor, direction);
end
for fieldIndex = 1:numel(muryFields)
    newGrid = copyGridMap(newGrid, oldGrid, muryFields{fieldIndex}, ...
        'mury', oldSpec, newSpec, factor, direction);
end
end

function newGrid = copyGridMap(newGrid, oldGrid, fieldName, layout, ...
        oldSpec, newSpec, factor, direction)
if ~isfield(oldGrid, fieldName)
    return
end
oldMap = oldGrid.(fieldName);
validateMap(oldMap, layoutSize(oldSpec, layout), ['cfg.grid.' fieldName]);
newGrid.(fieldName) = resampleMap( ...
    oldMap, layout, oldSpec, newSpec, factor, direction, 'nearest');
end

function newPml = resizePml(oldPml, oldSpec, newSpec, factor, direction)
if ~isstruct(oldPml) || ~isscalar(oldPml)
    error('fdtdmat:resampleCfg:InvalidPml', 'cfg.pml must be a scalar struct.');
end
newPml = oldPml;
fieldNames = fieldnames(oldPml);
for fieldIndex = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIndex};
    value = oldPml.(fieldName);
    if (isnumeric(value) || islogical(value)) && ...
            isequal(size(value), [oldSpec.Nx oldSpec.Ny])
        if islogical(value)
            newPml.(fieldName) = logical(resampleMap( ...
                value, 'cell', oldSpec, newSpec, factor, direction, ...
                'nearest', false));
        else
            newPml.(fieldName) = resampleMap( ...
                value, 'cell', oldSpec, newSpec, factor, direction, ...
                'linear', false);
        end
    end
end

indexScale = oldSpec.dx / newSpec.dx;
usesPhysicalUnits = isfield(oldPml, 'units') && ...
    strcmpi(char(string(oldPml.units)), 'physical');
if isfield(oldPml, 'thickness') && isnumeric(oldPml.thickness) && ...
        isscalar(oldPml.thickness)
    if usesPhysicalUnits
        newPml.thickness = oldPml.thickness;
    else
        newPml.thickness = max(1, round(oldPml.thickness * indexScale));
    end
end
if ~usesPhysicalUnits
    if isfield(oldPml, 'center') && isnumeric(oldPml.center) && ...
            numel(oldPml.center) == 2
        newPml.center = mapIndexCoordinates( ...
            oldPml.center, oldSpec, newSpec, false);
    end
    scalarLengths = {'outerRadius', 'innerRadius'};
    for lengthIndex = 1:numel(scalarLengths)
        fieldName = scalarLengths{lengthIndex};
        if isfield(oldPml, fieldName) && isnumeric(oldPml.(fieldName)) && ...
                isscalar(oldPml.(fieldName))
            newPml.(fieldName) = oldPml.(fieldName) * indexScale;
        end
    end
    depthFields = {'depthX', 'depthY', 'depth', 'radialDistance'};
    for depthIndex = 1:numel(depthFields)
        fieldName = depthFields{depthIndex};
        if isfield(newPml, fieldName) && isnumeric(newPml.(fieldName)) && ...
                isequal(size(newPml.(fieldName)), [newSpec.Nx newSpec.Ny])
            newPml.(fieldName) = newPml.(fieldName) * indexScale;
        end
    end
end
end

function newAntennas = resizeAntennas( ...
        oldAntennas, oldSpec, newSpec, factor, direction)
if ~isstruct(oldAntennas) || ~isscalar(oldAntennas)
    error('fdtdmat:resampleCfg:InvalidAntennas', ...
        'cfg.antennas must be a scalar struct.');
end
newAntennas = oldAntennas;
if isfield(oldAntennas, 'pos') && ~isempty(oldAntennas.pos)
    positions = mapIndexCoordinates( ...
        oldAntennas.pos, oldSpec, newSpec, true);
    validatePositions(positions, newSpec);
    if size(unique(positions, 'rows'), 1) ~= size(positions, 1)
        error('fdtdmat:resampleCfg:AntennaCollision', ...
            'Resampling maps multiple antennas to the same grid cell.');
    end
    newAntennas.pos = positions;
end
if isfield(oldAntennas, 'center') && numel(oldAntennas.center) == 2
    newAntennas.center = mapIndexCoordinates( ...
        oldAntennas.center, oldSpec, newSpec, false);
end
lengthFields = {'radius', 'pmlPadding', 'focusPadding'};
for fieldIndex = 1:numel(lengthFields)
    fieldName = lengthFields{fieldIndex};
    if isfield(oldAntennas, fieldName) && isnumeric(oldAntennas.(fieldName)) && ...
            isscalar(oldAntennas.(fieldName))
        scaledValue = round(oldAntennas.(fieldName) * oldSpec.dx / newSpec.dx);
        if strcmp(fieldName, 'radius')
            scaledValue = max(1, scaledValue);
        else
            scaledValue = max(0, scaledValue);
        end
        newAntennas.(fieldName) = scaledValue;
    end
end
boundFields = {'xMin', 'xMax', 'yMin', 'yMax'};
for fieldIndex = 1:numel(boundFields)
    fieldName = boundFields{fieldIndex};
    if ~isfield(oldAntennas, fieldName) || ...
            ~isnumeric(oldAntennas.(fieldName)) || ...
            ~isscalar(oldAntennas.(fieldName))
        continue
    end
    if startsWith(fieldName, 'x')
        oldSpacing = oldSpec.dx;
        newSpacing = newSpec.dx;
    else
        oldSpacing = oldSpec.dy;
        newSpacing = newSpec.dy;
    end
    newAntennas.(fieldName) = round( ...
        (oldAntennas.(fieldName) - 1) * oldSpacing / newSpacing + 1);
end
if isfield(oldAntennas, 'doiMask') && ~isempty(oldAntennas.doiMask)
    validateMap(oldAntennas.doiMask, [oldSpec.Nx oldSpec.Ny], ...
        'cfg.antennas.doiMask');
    newAntennas.doiMask = logical(resampleMap( ...
        oldAntennas.doiMask, 'cell', oldSpec, newSpec, ...
        factor, direction, 'nearest'));
    if ~any(newAntennas.doiMask(:))
        error('fdtdmat:resampleCfg:EmptyDoi', ...
            'Resampling produced an empty antenna DOI mask.');
    end
end
end

function newTargets = resizeTargets( ...
        oldTargets, oldSpec, newSpec, factor, direction)
if ~isstruct(oldTargets)
    error('fdtdmat:resampleCfg:InvalidTargets', ...
        'cfg.targets must be a struct array.');
end
newTargets = oldTargets;
for targetIndex = 1:numel(oldTargets)
    oldTarget = oldTargets(targetIndex);
    if isfield(oldTarget, 'mask') && ~isempty(oldTarget.mask)
        validateMap(oldTarget.mask, [oldSpec.Nx oldSpec.Ny], ...
            sprintf('cfg.targets(%d).mask', targetIndex));
        newTargets(targetIndex).mask = logical(resampleMap( ...
            oldTarget.mask, 'cell', oldSpec, newSpec, ...
            factor, direction, 'nearest'));
    end
    if ~isfield(oldTarget, 'properties') || ...
            ~isstruct(oldTarget.properties) || ...
            ~isfield(oldTarget, 'name')
        continue
    end
    properties = oldTarget.properties;
    switch lower(char(string(oldTarget.name)))
        case 'circle'
            if isfield(properties, 'center')
                properties.center = mapIndexCoordinates( ...
                    properties.center, oldSpec, newSpec, false);
            end
            if isfield(properties, 'radius')
                properties.radius = properties.radius * oldSpec.dx / newSpec.dx;
            end
        case 'rectangle'
            if isfield(properties, 'bounds') && numel(properties.bounds) == 4
                bounds = double(properties.bounds(:).');
                xBounds = (bounds(1:2) - 1) * oldSpec.dx / newSpec.dx + 1;
                yBounds = (bounds(3:4) - 1) * oldSpec.dy / newSpec.dy + 1;
                properties.bounds = [xBounds yBounds];
            end
        case 'triangle'
            if isfield(properties, 'vertices')
                properties.vertices = mapIndexCoordinates( ...
                    properties.vertices, oldSpec, newSpec, false);
            end
    end
    newTargets(targetIndex).properties = properties;
end
end

function outputMap = resampleMap(inputMap, layout, oldSpec, newSpec, ...
        factor, direction, upMethod, filterOnDownsample)
if nargin < 8
    filterOnDownsample = true;
end
[sourceX, sourceY] = layoutCoordinates(oldSpec, layout);
[targetX, targetY] = layoutCoordinates(newSpec, layout);
workingMap = double(inputMap);
if strcmp(direction, 'down')
    if filterOnDownsample
        workingMap = boxSmoothReplicate(workingMap, factor);
        method = 'linear';
    else
        method = upMethod;
    end
else
    method = upMethod;
end
targetX = min(max(targetX, sourceX(1)), sourceX(end));
targetY = min(max(targetY, sourceY(1)), sourceY(end));
if isscalar(sourceX) || isscalar(sourceY)
    outputMap = separableInterpolate( ...
        workingMap, sourceX, sourceY, targetX, targetY, method);
else
    [targetYGrid, targetXGrid] = meshgrid(targetY, targetX);
    outputMap = interp2( ...
        sourceY, sourceX, workingMap, targetYGrid, targetXGrid, method);
end
if islogical(inputMap)
    outputMap = outputMap >= 0.5;
end
end

function outputMap = separableInterpolate( ...
        inputMap, sourceX, sourceY, targetX, targetY, method)
if isscalar(sourceX)
    rowResampled = repmat(inputMap(1, :), numel(targetX), 1);
else
    rowResampled = interp1(sourceX, inputMap, targetX, method);
end
if isscalar(sourceY)
    outputMap = repmat(rowResampled(:, 1), 1, numel(targetY));
else
    outputMap = interp1(sourceY, rowResampled.', targetY, method).';
end
end

function smoothedMap = boxSmoothReplicate(inputMap, factor)
padBefore = floor((factor - 1) / 2);
padAfter = factor - 1 - padBefore;
rowIndex = [ones(1, padBefore), 1:size(inputMap, 1), ...
    size(inputMap, 1) * ones(1, padAfter)];
columnIndex = [ones(1, padBefore), 1:size(inputMap, 2), ...
    size(inputMap, 2) * ones(1, padAfter)];
paddedMap = inputMap(rowIndex, columnIndex);
smoothedMap = conv2(paddedMap, ones(factor) / factor^2, 'valid');
end

function [xCoordinates, yCoordinates] = layoutCoordinates(spec, layout)
switch layout
    case 'cell'
        xCoordinates = (0:spec.Nx-1) .* spec.dx;
        yCoordinates = (0:spec.Ny-1) .* spec.dy;
    case 'murx'
        xCoordinates = (0:spec.Nx-1) .* spec.dx;
        yCoordinates = (0.5:spec.Ny-1.5) .* spec.dy;
    case 'mury'
        xCoordinates = (0.5:spec.Nx-1.5) .* spec.dx;
        yCoordinates = (0:spec.Ny-1) .* spec.dy;
    otherwise
        error('fdtdmat:resampleCfg:UnknownLayout', ...
            'Unknown grid-map layout "%s".', layout);
end
end

function mapSize = layoutSize(spec, layout)
switch layout
    case 'cell'
        mapSize = [spec.Nx spec.Ny];
    case 'murx'
        mapSize = [spec.Nx spec.Ny - 1];
    case 'mury'
        mapSize = [spec.Nx - 1 spec.Ny];
end
end

function mappedCoordinates = mapIndexCoordinates( ...
        coordinates, oldSpec, newSpec, roundResult)
if ~isnumeric(coordinates) || size(coordinates, 2) ~= 2 || ...
        any(~isfinite(coordinates(:)))
    error('fdtdmat:resampleCfg:InvalidCoordinates', ...
        'Index coordinates must be a finite numeric array with two columns.');
end
mappedCoordinates = (double(coordinates) - 1) .* ...
    [oldSpec.dx oldSpec.dy] ./ [newSpec.dx newSpec.dy] + 1;
if roundResult
    mappedCoordinates = round(mappedCoordinates);
end
end

function validatePositions(positions, spec)
if any(positions(:, 1) < 1 | positions(:, 1) > spec.Nx | ...
        positions(:, 2) < 1 | positions(:, 2) > spec.Ny)
    error('fdtdmat:resampleCfg:AntennaOutOfBounds', ...
        'Resampled antenna positions must remain inside the grid.');
end
end

function validateMap(value, expectedSize, qualifiedName)
if ~(isnumeric(value) || islogical(value)) || ~ismatrix(value) || ...
        ~isequal(size(value), expectedSize) || ...
        (isnumeric(value) && any(~isfinite(value(:))))
    error('fdtdmat:resampleCfg:InvalidMap', ...
        '%s must be a finite 2-D map of size [%d %d].', ...
        qualifiedName, expectedSize(1), expectedSize(2));
end
end

function value = positiveInteger(value, errorId)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 1 || value ~= round(value)
    error(errorId, 'Expected a positive integer scalar.');
end
value = double(value);
end

function value = positiveScalar(value, errorId)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    error(errorId, 'Expected a positive finite numeric scalar.');
end
value = double(value);
end

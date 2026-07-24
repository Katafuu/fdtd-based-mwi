function cfg = buildCircularAntennaArrayIdx(cfg)
%buildCircularAntennaArrayIdx Add a rounded circular antenna array to cfg.

validateCfgField(cfg, 'grid');
validateCfgField(cfg, 'antennas');

gridSize = rowVector(cfg.grid.sizeXY, 2, 'buildCircularAntennaArrayIdx:InvalidGridSize');
numAntennas = positiveInteger(cfg.antennas.numAntennas, ...
    'buildCircularAntennaArrayIdx:InvalidNumAntennas');
center = rowVector(cfg.antennas.center, 2, ...
    'buildCircularAntennaArrayIdx:InvalidCenter');
radius = positiveInteger(cfg.antennas.radius, ...
    'buildCircularAntennaArrayIdx:InvalidRadius');
focusPadding = nonnegativeInteger(cfg.antennas.focusPadding, ...
    'buildCircularAntennaArrayIdx:InvalidFocusPadding');

positions = zeros(numAntennas, 2);
for k = 1:numAntennas
    angle = 2*pi*(k - 1)/numAntennas;
    positions(k, 1) = center(1) + radius*cos(angle);
    positions(k, 2) = center(2) + radius*sin(angle);
    positions(k, :) = round(positions(k, :));
end

if any(positions(:, 1) < 1 | positions(:, 1) > gridSize(1) | ...
        positions(:, 2) < 1 | positions(:, 2) > gridSize(2))
    error('buildCircularAntennaArrayIdx:AntennaOutOfBounds', ...
        'Circular antenna positions must lie inside cfg.grid.sizeXY.');
end
if size(unique(positions, 'rows'), 1) ~= numAntennas
    error('buildCircularAntennaArrayIdx:DuplicateAntennaPositions', ...
        'Rounded circular antenna positions must be unique.');
end

doiRadius = radius - focusPadding;
if doiRadius <= 0
    error('buildCircularAntennaArrayIdx:InvalidDOIRadius', ...
        'cfg.antennas.radius must be larger than cfg.antennas.focusPadding.');
end

[focusX, focusY] = ndgrid(1:gridSize(1), 1:gridSize(2));
doiMask = hypot(focusX - center(1), focusY - center(2)) <= doiRadius;

cfg.antennas.numAntennas = numAntennas;
cfg.antennas.center = center;
cfg.antennas.radius = radius;
cfg.antennas.pos = positions;
cfg.antennas.doiMask = doiMask;
cfg.antennas.focusPadding = focusPadding;
cfg.antennas.arrayType = 'circular';
end

function validateCfgField(cfg, fieldName)
if ~isstruct(cfg) || ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    error('buildCircularAntennaArrayIdx:MissingCfgField', ...
        'cfg.%s is required.', fieldName);
end
end

function value = rowVector(value, expectedLength, errorId)
value = double(value(:).');
if numel(value) ~= expectedLength || any(~isfinite(value)) || ...
        any(value ~= round(value))
    error(errorId, 'Value must contain %d finite integer entries.', expectedLength);
end
end

function value = positiveInteger(value, errorId)
if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ...
        value ~= round(value) || value < 1
    error(errorId, 'Value must be a positive integer scalar.');
end
value = double(value);
end

function value = nonnegativeInteger(value, errorId)
if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ...
        value ~= round(value) || value < 0
    error(errorId, 'Value must be a nonnegative integer scalar.');
end
value = double(value);
end

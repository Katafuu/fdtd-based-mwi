function cfg = buildPlanarAntennaIdx(cfg)
%buildPlanarAntennaIdx Add a top-side planar antenna line to cfg.

validateCfgField(cfg, 'grid');
validateCfgField(cfg, 'pml');
validateCfgField(cfg, 'antennas');

gridSize = rowVector(cfg.grid.sizeXY, 2, 'buildPlanarAntennaIdx:InvalidGridSize');
numAntennas = positiveInteger(cfg.antennas.numAntennas, ...
    'buildPlanarAntennaIdx:InvalidNumAntennas');
pmlThickness = nonnegativeInteger(cfg.pml.thickness, ...
    'buildPlanarAntennaIdx:InvalidPmlThickness');
pmlPadding = nonnegativeInteger(cfg.antennas.pmlPadding, ...
    'buildPlanarAntennaIdx:InvalidPmlPadding');
focusPadding = nonnegativeInteger(cfg.antennas.focusPadding, ...
    'buildPlanarAntennaIdx:InvalidFocusPadding');

orientation = char(cfg.antennas.orientation);
if ~strcmpi(orientation, 'top')
    error('buildPlanarAntennaIdx:UnsupportedOrientation', ...
        'Only cfg.antennas.orientation = ''top'' is currently supported.');
end
orientation = 'top';

Nx = gridSize(1);
Ny = gridSize(2);
antennaY = 1 + pmlThickness + pmlPadding;
xMin = 1 + pmlThickness + pmlPadding;
xMax = Nx - pmlThickness - pmlPadding;

if xMin > xMax
    error('buildPlanarAntennaIdx:InvalidHorizontalExtent', ...
        'PML thickness and padding leave no horizontal room for antennas.');
end
if antennaY < 1 || antennaY > Ny || antennaY > Ny - pmlThickness
    error('buildPlanarAntennaIdx:InvalidAntennaY', ...
        'The planar antenna line must lie inside the non-PML grid.');
end

antennaX = round(linspace(xMin, xMax, numAntennas));
if numel(unique(antennaX)) ~= numAntennas
    error('buildPlanarAntennaIdx:DuplicateAntennaPositions', ...
        'Rounded planar antenna positions must be unique.');
end

positions = [antennaX(:), antennaY * ones(numAntennas, 1)];

doiXMin = xMin + focusPadding;
doiXMax = xMax - focusPadding;
doiYMin = antennaY + focusPadding;
doiYMax = Ny - pmlThickness - pmlPadding - focusPadding;
if doiXMin > doiXMax || doiYMin > doiYMax
    error('buildPlanarAntennaIdx:InvalidDOIExtent', ...
        'PML padding and focus padding leave no valid planar DOI.');
end

doiMask = false(Nx, Ny);
doiMask(doiXMin:doiXMax, doiYMin:doiYMax) = true;

center = [mean(antennaX), antennaY];
radius = max(abs(antennaX - center(1)));

cfg.antennas.numAntennas = numAntennas;
cfg.antennas.center = center;
cfg.antennas.radius = radius;
cfg.antennas.pos = positions;
cfg.antennas.doiMask = doiMask;
cfg.antennas.focusPadding = focusPadding;
cfg.antennas.arrayType = 'planar';
cfg.antennas.orientation = orientation;
cfg.antennas.pmlPadding = pmlPadding;
end

function validateCfgField(cfg, fieldName)
if ~isstruct(cfg) || ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    error('buildPlanarAntennaIdx:MissingCfgField', ...
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

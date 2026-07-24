function grid = createGrid(sizeXY, spacingXY, background, originPhysical)
%createGrid Create grid coordinates and default material maps.

if nargin < 2 || isempty(spacingXY)
    spacingXY = [1 1];
end
if nargin < 3 || isempty(background)
    background = struct();
end
if nargin < 4 || isempty(originPhysical)
    originPhysical = [0 0];
end

sizeXY = rowVector(sizeXY, 2, 'fdtdmat:createGrid:InvalidSize');
spacingXY = rowVector(spacingXY, 2, 'fdtdmat:createGrid:InvalidSpacing');
originPhysical = rowVector(originPhysical, 2, 'fdtdmat:createGrid:InvalidOrigin');
background = fillBackgroundDefaults(background);

if any(sizeXY < 1) || any(fix(sizeXY) ~= sizeXY)
    error('fdtdmat:createGrid:InvalidSize', 'sizeXY must contain positive integer grid dimensions.');
end
if any(spacingXY <= 0)
    error('fdtdmat:createGrid:InvalidSpacing', 'spacingXY must contain positive values.');
end

[xIndex, yIndex] = ndgrid(1:sizeXY(1), 1:sizeXY(2));

grid = struct();
grid.sizeXY = sizeXY;
grid.spacingXY = spacingXY;
grid.originPhysical = originPhysical;
grid.xIndex = xIndex;
grid.yIndex = yIndex;
grid.xPhysical = originPhysical(1) + (xIndex - 1) * spacingXY(1);
grid.yPhysical = originPhysical(2) + (yIndex - 1) * spacingXY(2);
grid.background = background;

grid.epsr = background.epsr * ones(sizeXY);
grid.murx = background.murx * ones(sizeXY(1), sizeXY(2) - 1);
grid.mury = background.mury * ones(sizeXY(1) - 1, sizeXY(2));
grid.cond_e = background.cond_e * ones(sizeXY);
grid.cond_m = background.cond_m * ones(sizeXY);
grid.epsr_bg = grid.epsr;
grid.murx_bg = grid.murx;
grid.mury_bg = grid.mury;
grid.cond_e_bg = grid.cond_e;
grid.cond_m_bg = grid.cond_m;
end

function value = rowVector(value, expectedLength, errorId)
value = double(value(:).');
if numel(value) ~= expectedLength || any(~isfinite(value))
    error(errorId, 'Expected a finite vector with %d elements.', expectedLength);
end
end

function background = fillBackgroundDefaults(background)
if ~isstruct(background)
    error('fdtdmat:createGrid:InvalidBackground', 'background must be a struct.');
end

background = withDefault(background, 'epsr', 1.0);
background = withDefault(background, 'murx', 1.0);
background = withDefault(background, 'mury', 1.0);
background = withDefault(background, 'cond_e', 0.0);
background = withDefault(background, 'cond_m', 0.0);

validateBackgroundScalar(background.epsr, 'epsr');
validateBackgroundScalar(background.murx, 'murx');
validateBackgroundScalar(background.mury, 'mury');
validateBackgroundScalar(background.cond_e, 'cond_e');
validateBackgroundScalar(background.cond_m, 'cond_m');

if background.epsr <= 0 || background.murx <= 0 || background.mury <= 0
    error('fdtdmat:createGrid:InvalidBackground', 'background epsr, murx, and mury must be positive.');
end
if background.cond_e < 0 || background.cond_m < 0
    error('fdtdmat:createGrid:InvalidBackground', 'background conductivities must be nonnegative.');
end
end

function background = withDefault(background, fieldName, defaultValue)
if ~isfield(background, fieldName) || isempty(background.(fieldName))
    background.(fieldName) = defaultValue;
end
end

function validateBackgroundScalar(value, fieldName)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    error('fdtdmat:createGrid:InvalidBackground', ...
        'background.%s must be a finite numeric scalar.', fieldName);
end
end

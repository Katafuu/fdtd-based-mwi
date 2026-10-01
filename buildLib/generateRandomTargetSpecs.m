function targets = generateRandomTargetSpecs(grid, doiMask, opts, occupiedMask)
%generateRandomTargetSpecs Generate random target geometry/material specs.
% occupiedMask optionally reserves cells for previously placed exact targets.

if nargin < 3 || isempty(opts)
    opts = struct();
end

gridSize = double(grid.sizeXY(:).');
if numel(gridSize) ~= 2 || any(gridSize < 1) || any(gridSize ~= round(gridSize))
    error('generateRandomTargetSpecs:InvalidGridSize', ...
        'grid.sizeXY must contain two positive integer dimensions.');
end

doiMask = logical(doiMask);
if ~isequal(size(doiMask), gridSize)
    error('generateRandomTargetSpecs:InvalidDOIMask', ...
        'doiMask must match grid.sizeXY.');
end
if ~any(doiMask(:))
    error('generateRandomTargetSpecs:EmptyDOIMask', ...
        'doiMask must contain at least one true pixel.');
end

numTargetsRange = integerRange(getOption(opts, 'numTargetsRange', [1 3]), ...
    'generateRandomTargetSpecs:InvalidNumTargetsRange');
allowedShapes = string(getOption(opts, 'allowedShapes', ["circle" "rectangle" "triangle"]));
allowedShapes = lower(allowedShapes(:).');
if isempty(allowedShapes) || any(~ismember(allowedShapes, ["circle" "rectangle" "triangle"]))
    error('generateRandomTargetSpecs:InvalidAllowedShapes', ...
        'opts.allowedShapes may contain only "circle", "rectangle", and "triangle".');
end

radiusRange = integerRange(getOption(opts, 'radiusRange', [10 30]), ...
    'generateRandomTargetSpecs:InvalidRadiusRange');
sideRange = integerRange(getOption(opts, 'sideRange', [20 50]), ...
    'generateRandomTargetSpecs:InvalidSideRange');
epsrRange = numericRange(getOption(opts, 'epsrRange', [2 6]), ...
    'generateRandomTargetSpecs:InvalidEpsrRange');
condRange = numericRange(getOption(opts, 'condRange', [0.02 0.15]), ...
    'generateRandomTargetSpecs:InvalidCondRange');
maxAttempts = positiveInteger(getOption(opts, 'maxAttempts', 500), ...
    'generateRandomTargetSpecs:InvalidMaxAttempts');
allowOverlap = logicalScalar(getOption(opts, 'allowOverlap', false), ...
    'generateRandomTargetSpecs:InvalidAllowOverlap');

numTargets = randi(numTargetsRange);
[candidateX, candidateY] = find(doiMask);
if nargin < 4 || isempty(occupiedMask)
    occupiedMask = false(gridSize);
else
    if ~isequal(size(occupiedMask), gridSize)
        error('generateRandomTargetSpecs:InvalidOccupiedMask', ...
            'occupiedMask must match grid.sizeXY.');
    end
    occupiedMask = logical(occupiedMask);
end
targetTemplate = struct( ...
    'mask', [], ...
    'name', "", ...
    'material', struct('epsr', NaN, 'murx', NaN, 'mury', NaN, 'cond_e', NaN, 'cond_m', NaN), ...
    'properties', struct());
targets = repmat(targetTemplate, 1, numTargets);

for targetIdx = 1:numTargets
    placedTarget = false;
    for attemptIdx = 1:maxAttempts
        name = allowedShapes(randi(numel(allowedShapes)));
        candidateIdx = randi(numel(candidateX));
        center = [candidateX(candidateIdx), candidateY(candidateIdx)];
        target = buildCandidateTarget(targetTemplate, name, center, ...
            radiusRange, sideRange, epsrRange, condRange);
        region = targetRegion(grid, target);
        targetMask = region.mask;

        if ~any(targetMask(:))
            continue
        end
        if any(targetMask(:) & ~doiMask(:))
            continue
        end
        if ~allowOverlap && any(targetMask(:) & occupiedMask(:))
            continue
        end

        target.mask = targetMask;
        targets(targetIdx) = target;
        occupiedMask = occupiedMask | targetMask;
        placedTarget = true;
        break
    end

    if ~placedTarget
        error('generateRandomTargetSpecs:PlacementFailed', ...
            'Could not place target %d inside the DOI after %d attempts.', ...
            targetIdx, maxAttempts);
    end
end
end

function target = buildCandidateTarget(template, name, center, radiusRange, sideRange, epsrRange, condRange)
target = template;
target.name = name;
target.material = struct( ...
    'epsr', uniformRange(epsrRange), ...
    'murx', 1.0, ...
    'mury', 1.0, ...
    'cond_e', uniformRange(condRange), ...
    'cond_m', 0.0);

switch name
    case "circle"
        target.properties = struct( ...
            'center', double(center), ...
            'radius', randi(radiusRange));
    case "rectangle"
        side = randi(sideRange);
        leftExtent = floor((side - 1) / 2);
        rightExtent = ceil((side - 1) / 2);
        target.properties = struct('bounds', [ ...
            center(1) - leftExtent, center(1) + rightExtent, ...
            center(2) - leftExtent, center(2) + rightExtent]);
    case "triangle"
        side = randi(sideRange);
        halfSide = 0.5 * side;
        height = sqrt(3) * side / 2;
        target.properties = struct('vertices', [ ...
            center(1), center(2) - 2 * height / 3; ...
            center(1) - halfSide, center(2) + height / 3; ...
            center(1) + halfSide, center(2) + height / 3]);
end
end

function region = targetRegion(grid, target)
switch target.name
    case "circle"
        region = fdtdgeom.shape_circle(grid, target.properties.center, ...
            target.properties.radius, 'index');
    case "rectangle"
        region = fdtdgeom.shape_rectangle(grid, target.properties.bounds, 'index');
    case "triangle"
        region = fdtdgeom.shape_polygon(grid, target.properties.vertices, 'index');
    otherwise
        error('generateRandomTargetSpecs:UnsupportedShape', ...
            'Unsupported target name "%s".', target.name);
end
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end

function range = integerRange(value, errorId)
range = double(value(:).');
if numel(range) ~= 2 || any(~isfinite(range)) || any(range ~= round(range)) || ...
        range(1) < 1 || range(2) < range(1)
    error(errorId, 'Expected a two-element positive integer range [min max].');
end
end

function range = numericRange(value, errorId)
range = double(value(:).');
if numel(range) ~= 2 || any(~isfinite(range)) || range(2) < range(1)
    error(errorId, 'Expected a two-element numeric range [min max].');
end
end

function value = positiveInteger(value, errorId)
if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ...
        value ~= round(value) || value < 1
    error(errorId, 'Expected a positive integer scalar.');
end
value = double(value);
end

function value = logicalScalar(value, errorId)
if ~isscalar(value) || ~(islogical(value) || isnumeric(value))
    error(errorId, 'Expected a logical scalar.');
end
value = logical(value);
end

function value = uniformRange(range)
value = range(1) + rand() * (range(2) - range(1));
end
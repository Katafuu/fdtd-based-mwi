function [solution, info] = solve_tsvd(operator, rhs, options)
%solve_tsvd Solve a linear system with an energy-truncated SVD.

if nargin < 3 || isempty(options)
    options = struct();
end
if ~isnumeric(operator) || ~ismatrix(operator) || isempty(operator) || ...
        any(~isfinite(operator(:)))
    error('solve_tsvd:InvalidOperator', ...
        'operator must be a nonempty finite numeric matrix.');
end
if ~isnumeric(rhs) || ~isvector(rhs) || size(operator, 1) ~= numel(rhs) || ...
        any(~isfinite(rhs(:)))
    error('solve_tsvd:InvalidRightHandSide', ...
        'rhs must be a finite vector matching the operator row count.');
end
rhs = cast(rhs(:), 'like', operator);

energyThreshold = getScalarOption(options, 'energyThreshold', 0.99, 0, 1);
initialRank = getIntegerOption(options, 'initialRank', 32);
maximumRank = getIntegerOption(options, 'maximumRank', 512);
fullSvdMaxDimension = getIntegerOption(options, 'fullSvdMaxDimension', 512);
relativeSingularTolerance = getScalarOption(options, ...
    'relativeSingularTolerance', 10*eps(class(operator)), 0, 1);

minimumDimension = min(size(operator));
maximumRank = min(maximumRank, minimumDimension);
totalEnergy = double(norm(operator, 'fro'))^2;
if totalEnergy == 0
    solution = zeros(size(operator, 2), 1);
    info = emptyInfo(energyThreshold, totalEnergy);
    return
end

usedFullSvd = minimumDimension <= fullSvdMaxDimension;
if usedFullSvd
    [leftVectors, singularMatrix, rightVectors] = svd(operator, 'econ');
    singularValues = diag(singularMatrix);
    candidateCount = min(maximumRank, numel(singularValues));
    leftVectors = leftVectors(:, 1:candidateCount);
    rightVectors = rightVectors(:, 1:candidateCount);
    singularValues = singularValues(1:candidateCount);
else
    maximumIterativeRank = min(maximumRank, minimumDimension - 1);
    if maximumIterativeRank < 1
        error('solve_tsvd:InvalidRank', ...
            'The iterative TSVD path requires a minimum dimension above one.');
    end
    trialRank = min(initialRank, maximumIterativeRank);
    while true
        [leftVectors, singularMatrix, rightVectors] = ...
            svds(operator, trialRank);
        singularValues = diag(singularMatrix);
        [singularValues, order] = sort(singularValues, 'descend');
        leftVectors = leftVectors(:, order);
        rightVectors = rightVectors(:, order);

        trialEnergy = sum(double(abs(singularValues)).^2) / totalEnergy;
        if trialEnergy >= energyThreshold || trialRank >= maximumIterativeRank
            break
        end
        trialRank = min(maximumIterativeRank, max(trialRank + 1, 2*trialRank));
    end
end

cumulativeEnergy = cumsum(double(abs(singularValues)).^2) / totalEnergy;
retainedRank = find(cumulativeEnergy >= energyThreshold, 1, 'first');
if isempty(retainedRank)
    retainedRank = numel(singularValues);
end

relativeMagnitude = abs(singularValues(1:retainedRank)) / ...
    max(abs(singularValues(1)), realmin(class(singularValues)));
numericalRank = find(relativeMagnitude > relativeSingularTolerance, 1, 'last');
if isempty(numericalRank)
    error('solve_tsvd:NumericallyZeroOperator', ...
        'No singular value exceeds the configured numerical tolerance.');
end
retainedRank = numericalRank;

retainedSingularValues = singularValues(1:retainedRank);
leftVectors = leftVectors(:, 1:retainedRank);
rightVectors = rightVectors(:, 1:retainedRank);
coefficients = (leftVectors' * rhs) ./ retainedSingularValues;
solution = double(rightVectors * coefficients);

achievedEnergy = sum(double(abs(retainedSingularValues)).^2) / totalEnergy;
info = struct();
info.retainedRank = retainedRank;
info.achievedEnergy = achievedEnergy;
info.requestedEnergy = energyThreshold;
info.totalEnergy = totalEnergy;
info.singularValues = double(retainedSingularValues);
info.conditionEstimate = double(abs(retainedSingularValues(1)) / ...
    abs(retainedSingularValues(end)));
info.hitRankCap = achievedEnergy < energyThreshold;
info.usedFullSvd = usedFullSvd;
end

function info = emptyInfo(energyThreshold, totalEnergy)
info = struct( ...
    'retainedRank', 0, ...
    'achievedEnergy', 0, ...
    'requestedEnergy', energyThreshold, ...
    'totalEnergy', totalEnergy, ...
    'singularValues', zeros(0, 1), ...
    'conditionEstimate', NaN, ...
    'hitRankCap', false, ...
    'usedFullSvd', true);
end

function value = getScalarOption(options, name, defaultValue, lowerBound, upperBound)
value = getOption(options, name, defaultValue);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value <= lowerBound || value > upperBound
    error('solve_tsvd:InvalidOption', ...
        'options.%s must lie in (%g, %g].', name, lowerBound, upperBound);
end
value = double(value);
end

function value = getIntegerOption(options, name, defaultValue)
value = getOption(options, name, defaultValue);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 1 || value ~= round(value)
    error('solve_tsvd:InvalidOption', ...
        'options.%s must be a positive integer.', name);
end
value = double(value);
end

function value = getOption(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

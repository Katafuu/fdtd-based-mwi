function [operator, rhs, metadata] = buildDbimLinearSystem( ...
    measuredFields, predictedFields, doiFields, pairs, analyticGreen, ...
    frequencies, cfg, options)
%buildDbimLinearSystem Assemble the real fixed-Green multi-frequency system.

if nargin < 8 || isempty(options)
    options = struct();
end
frequencies = double(frequencies(:).');
numFrequencies = numel(frequencies);
if isempty(frequencies) || any(~isfinite(frequencies)) || any(frequencies <= 0)
    error('buildDbimLinearSystem:InvalidFrequencies', ...
        'frequencies must contain finite positive values.');
end
if ~isnumeric(pairs) || size(pairs, 2) ~= 2 || isempty(pairs) || ...
        any(~isfinite(pairs(:))) || any(pairs(:) ~= round(pairs(:)))
    error('buildDbimLinearSystem:InvalidPairs', ...
        'pairs must be a nonempty integer N-by-2 array [tx rx].');
end

numAntennas = size(measuredFields, 1);
if size(measuredFields, 2) ~= numAntennas || ...
        size(measuredFields, 3) ~= numFrequencies || ...
        size(predictedFields, 1) ~= numAntennas || ...
        size(predictedFields, 2) ~= numAntennas || ...
        size(predictedFields, 3) ~= numFrequencies || ...
        any(~isfinite(measuredFields(:))) || any(~isfinite(predictedFields(:)))
    error('buildDbimLinearSystem:InvalidReceiverFields', ...
        'Receiver fields must be finite antenna-by-antenna-by-frequency arrays.');
end
numPixels = size(doiFields, 1);
if size(doiFields, 2) ~= numAntennas || ...
        size(doiFields, 3) ~= numFrequencies || ...
        any(~isfinite(doiFields(:)))
    error('buildDbimLinearSystem:InvalidDoiFields', ...
        'doiFields must be a finite pixel-by-antenna-by-frequency array.');
end
if ~isnumeric(analyticGreen) || ...
        size(analyticGreen, 1) ~= numPixels || ...
        size(analyticGreen, 2) ~= numAntennas || ...
        size(analyticGreen, 3) ~= numFrequencies || ...
        any(~isfinite(analyticGreen(:)))
    error('buildDbimLinearSystem:InvalidAnalyticGreen', ...
        ['analyticGreen must be a finite ' ...
         'pixel-by-antenna-by-frequency array.']);
end
if any(pairs(:) < 1) || any(pairs(:) > numAntennas) || ...
        any(pairs(:, 1) == pairs(:, 2))
    error('buildDbimLinearSystem:PairOutOfRange', ...
        'Pair indices must be distinct valid antenna indices.');
end
if ~isstruct(cfg) || ~all(isfield(cfg, {'mu0', 'eps0', 'dx', 'dy'}))
    error('buildDbimLinearSystem:InvalidConfig', ...
        'cfg must contain mu0, eps0, dx, and dy.');
end

epsrScale = positiveScalarOption(options, 'epsrScale', 1.0);
sigmaScale = positiveScalarOption(options, 'sigmaScale', 0.1);
blockSize = positiveIntegerOption(options, 'assemblyBlockSize', 64);
operatorPrecision = stringOption(options, 'operatorPrecision', "single");
if ~ismember(operatorPrecision, ["single" "double"])
    error('buildDbimLinearSystem:InvalidOperatorPrecision', ...
        'operatorPrecision must be "single" or "double".');
end

numPairs = size(pairs, 1);
numRows = 2 * numPairs * numFrequencies;
numColumns = 2 * numPixels;
operator = zeros(numRows, numColumns, char(operatorPrecision));
rhs = zeros(numRows, 1, char(operatorPrecision));
complexResidual = complex(zeros(numPairs, numFrequencies));
measuredPairFields = complex(zeros(numPairs, numFrequencies));
pairLinearIndex = sub2ind([numAntennas numAntennas], ...
    pairs(:, 2), pairs(:, 1));

for frequencyIndex = 1:numFrequencies
    omega = 2*pi*frequencies(frequencyIndex);
    measuredAtFrequency = measuredFields(:, :, frequencyIndex);
    predictedAtFrequency = predictedFields(:, :, frequencyIndex);
    measuredPairFields(:, frequencyIndex) = ...
        measuredAtFrequency(pairLinearIndex);
    complexResidual(:, frequencyIndex) = ...
        measuredAtFrequency(pairLinearIndex) - ...
        predictedAtFrequency(pairLinearIndex);

    rowOffset = (frequencyIndex - 1) * 2 * numPairs;
    realRows = rowOffset + (1:numPairs);
    imaginaryRows = rowOffset + numPairs + (1:numPairs);
    rhs(realRows) = cast(real(complexResidual(:, frequencyIndex)), ...
        char(operatorPrecision));
    rhs(imaginaryRows) = cast(imag(complexResidual(:, frequencyIndex)), ...
        char(operatorPrecision));

    physicalFactor = omega^2 * cfg.mu0 * cfg.eps0 * cfg.dx * cfg.dy;
    conductivityFactor = omega * cfg.eps0;
    for firstPair = 1:blockSize:numPairs
        pairBlock = firstPair:min(firstPair + blockSize - 1, numPairs);
        tx = pairs(pairBlock, 1);
        rx = pairs(pairBlock, 2);

        transmitterField = doiFields(:, tx, frequencyIndex).';
        receiverGreen = analyticGreen(:, rx, frequencyIndex).';
        complexOperatorBlock = physicalFactor .* ...
            (receiverGreen .* transmitterField);
        realBlock = real(complexOperatorBlock);
        imaginaryBlock = imag(complexOperatorBlock);

        operator(realRows(pairBlock), 1:numPixels) = cast( ...
            epsrScale .* realBlock, char(operatorPrecision));
        operator(imaginaryRows(pairBlock), 1:numPixels) = cast( ...
            epsrScale .* imaginaryBlock, char(operatorPrecision));
        operator(realRows(pairBlock), numPixels + (1:numPixels)) = cast( ...
            sigmaScale .* imaginaryBlock ./ conductivityFactor, ...
            char(operatorPrecision));
        operator(imaginaryRows(pairBlock), numPixels + (1:numPixels)) = cast( ...
            -sigmaScale .* realBlock ./ conductivityFactor, ...
            char(operatorPrecision));
    end
end

metadata = struct();
metadata.numPixels = numPixels;
metadata.numPairs = numPairs;
metadata.epsrScale = epsrScale;
metadata.sigmaScale = sigmaScale;
metadata.complexResidual = complexResidual;
metadata.measuredPairFields = measuredPairFields;
metadata.residualNorm = norm(complexResidual(:));
metadata.measuredNorm = norm(measuredPairFields(:));
metadata.normalizedResidual = metadata.residualNorm / ...
    max(metadata.measuredNorm, eps);
metadata.pairLinearIndex = pairLinearIndex;
metadata.greenFunction = "analytical homogeneous 2-D";
end

function value = positiveScalarOption(options, name, defaultValue)
value = getOption(options, name, defaultValue);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    error('buildDbimLinearSystem:InvalidOption', ...
        'options.%s must be a finite positive scalar.', name);
end
value = double(value);
end

function value = positiveIntegerOption(options, name, defaultValue)
value = positiveScalarOption(options, name, defaultValue);
if value ~= round(value)
    error('buildDbimLinearSystem:InvalidOption', ...
        'options.%s must be a positive integer.', name);
end
end

function value = stringOption(options, name, defaultValue)
value = string(getOption(options, name, defaultValue));
if ~isscalar(value)
    error('buildDbimLinearSystem:InvalidOption', ...
        'options.%s must be a string scalar.', name);
end
end

function value = getOption(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
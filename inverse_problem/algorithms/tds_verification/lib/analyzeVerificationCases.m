function results = analyzeVerificationCases(measurements)
%analyzeVerificationCases Compare baseline and maximum-chord estimates.

if isstruct(measurements)
    measurements = num2cell(measurements);
end
if ~iscell(measurements) || isempty(measurements)
    error('analyzeVerificationCases:InvalidMeasurements', ...
        'measurements must be a nonempty cell array or struct array.');
end

numCases = numel(measurements);
seed = nan(numCases, 1);
knownEpsr = nan(numCases, 1);
baselineEpsr = nan(numCases, 1);
baselineErrorPercent = nan(numCases, 1);
maximumChordEpsr = nan(numCases, 1);
maximumChordErrorPercent = nan(numCases, 1);
selectedDirectedRays = zeros(numCases, 1);
maximumChordMillimeters = nan(numCases, 1);

for caseIndex = 1:numCases
    measurement = measurements{caseIndex};
    if isfield(measurement, 'seed')
        seed(caseIndex) = measurement.seed;
    end
    knownEpsr(caseIndex) = measurement.targetEpsr;

    prominentDelaySteps = max(measurement.prominentDelaySteps, 0);
    deltaTime = prominentDelaySteps .* measurement.dt;
    averageEpsr = (sqrt(measurement.backgroundEpsr) + ...
        measurement.c0 .* deltaTime ./ measurement.distance).^2;
    [baselineRayEstimates, ~, ~, baselineAccepted] = ...
        evaluateTargetRays(averageEpsr, measurement.backgroundEpsr, ...
        measurement.antennaPositions, logical(measurement.targetMask), ...
        15, logical(measurement.pairEligible));
    baselineEpsr(caseIndex) = ...
        mean(baselineRayEstimates(baselineAccepted));

    [maximumChordEpsr(caseIndex), details] = ...
        estimateMaximumChordEpsr(measurement);
    selectedDirectedRays(caseIndex) = details.numSelectedDirectedRays;
    maximumChordMillimeters(caseIndex) = ...
        1e3 * details.maximumChordLength;

    baselineErrorPercent(caseIndex) = 100 * ...
        abs(baselineEpsr(caseIndex)-knownEpsr(caseIndex)) / ...
        knownEpsr(caseIndex);
    maximumChordErrorPercent(caseIndex) = 100 * ...
        abs(maximumChordEpsr(caseIndex)-knownEpsr(caseIndex)) / ...
        knownEpsr(caseIndex);
end

results = table(seed, knownEpsr, baselineEpsr, baselineErrorPercent, ...
    maximumChordEpsr, maximumChordErrorPercent, ...
    selectedDirectedRays, maximumChordMillimeters);
end

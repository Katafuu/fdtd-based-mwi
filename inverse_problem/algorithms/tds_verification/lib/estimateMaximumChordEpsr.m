function [recoveredEpsr, details] = estimateMaximumChordEpsr( ...
        measurement, options)
%estimateMaximumChordEpsr Estimate epsr from matched delays on longest chords.

if nargin < 2 || isempty(options)
    options = struct();
end
requiredFields = {'targetMask', 'antennaPositions', 'dx', 'dy', ...
    'pairEligible', 'matchedDelaySteps', 'c0', 'dt', 'backgroundEpsr'};
missingFields = requiredFields(~isfield(measurement, requiredFields));
if ~isempty(missingFields)
    error('estimateMaximumChordEpsr:MissingMeasurementField', ...
        'measurement is missing field(s): %s.', strjoin(missingFields, ', '));
end

relativeChordThreshold = getOption(options, ...
    'RelativeChordThreshold', 0.999);
sampleDensity = getOption(options, 'SampleDensity', 20);
validateattributes(relativeChordThreshold, {'numeric'}, ...
    {'real', 'finite', 'scalar', '>', 0, '<=', 1}, ...
    mfilename, 'RelativeChordThreshold');

chordLengths = computeCenterlineIntersectionLengths( ...
    logical(measurement.targetMask), measurement.antennaPositions, ...
    [measurement.dx measurement.dy], ...
    logical(measurement.pairEligible), sampleDensity);
delaySteps = measurement.matchedDelaySteps;
candidatePairs = logical(measurement.pairEligible) & ...
    isfinite(delaySteps) & delaySteps >= 0 & chordLengths > 0;
if ~any(candidatePairs, 'all')
    error('estimateMaximumChordEpsr:NoCandidateRays', ...
        'No finite matched-delay ray intersects the target centerline.');
end

maximumChordLength = max(chordLengths(candidatePairs));
selectedPairs = candidatePairs & ...
    chordLengths >= relativeChordThreshold*maximumChordLength;
opticalPathDifferences = measurement.c0 * measurement.dt .* delaySteps;
selectedLengths = chordLengths(selectedPairs);
selectedDifferences = opticalPathDifferences(selectedPairs);

refractiveIndexContrast = ...
    sum(selectedLengths .* selectedDifferences) / ...
    sum(selectedLengths.^2);
backgroundIndex = sqrt(measurement.backgroundEpsr);
recoveredIndex = backgroundIndex + refractiveIndexContrast;
recoveredEpsr = recoveredIndex^2;

perRayIndex = backgroundIndex + ...
    opticalPathDifferences ./ chordLengths;
perRayEpsr = perRayIndex.^2;
perRayEpsr(~candidatePairs) = NaN;

details = struct();
details.relativeChordThreshold = relativeChordThreshold;
details.sampleDensity = sampleDensity;
details.chordLengths = chordLengths;
details.maximumChordLength = maximumChordLength;
details.candidatePairs = candidatePairs;
details.selectedPairs = selectedPairs;
details.numSelectedDirectedRays = nnz(selectedPairs);
details.refractiveIndexContrast = refractiveIndexContrast;
details.recoveredIndex = recoveredIndex;
details.perRayEpsr = perRayEpsr;
details.selectedPerRayEpsr = perRayEpsr(selectedPairs);
end

function value = getOption(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

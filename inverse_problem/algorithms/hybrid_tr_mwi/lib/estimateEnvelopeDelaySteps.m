function [delaySteps, maximumScore] = estimateEnvelopeDelaySteps( ...
        referenceSignal, objectSignal, referenceArrivalStep, ...
        windowHalfWidthSteps, maximumDelaySteps)
%estimateEnvelopeDelaySteps Match a reference envelope at nonnegative lags.

validateattributes(referenceSignal, {'numeric'}, ...
    {'real', 'finite', 'vector', 'nonempty'}, mfilename, 'referenceSignal');
validateattributes(objectSignal, {'numeric'}, ...
    {'real', 'finite', 'vector', 'numel', numel(referenceSignal)}, ...
    mfilename, 'objectSignal');
validateattributes(referenceArrivalStep, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, ...
    mfilename, 'referenceArrivalStep');
validateattributes(windowHalfWidthSteps, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
    mfilename, 'windowHalfWidthSteps');
validateattributes(maximumDelaySteps, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'nonnegative'}, ...
    mfilename, 'maximumDelaySteps');

referenceEnvelope = abs(hilbert(referenceSignal(:)));
objectEnvelope = abs(hilbert(objectSignal(:)));
numSteps = numel(referenceEnvelope);
referenceCenter = round(referenceArrivalStep);
referenceIndices = max(1, referenceCenter-windowHalfWidthSteps) : ...
    min(numSteps, referenceCenter+windowHalfWidthSteps);
maximumUsableDelay = min(maximumDelaySteps, ...
    numSteps-referenceIndices(end));

if maximumUsableDelay < 0
    delaySteps = NaN;
    maximumScore = NaN;
    return
end

referenceSegment = referenceEnvelope(referenceIndices);
referenceNorm = norm(referenceSegment);
if referenceNorm <= eps(max(referenceEnvelope))
    delaySteps = NaN;
    maximumScore = NaN;
    return
end

scores = nan(maximumUsableDelay+1, 1);
for integerDelay = 0:maximumUsableDelay
    objectSegment = objectEnvelope(referenceIndices + integerDelay);
    objectNorm = norm(objectSegment);
    if objectNorm > eps(max(objectEnvelope))
        scores(integerDelay+1) = ...
            dot(referenceSegment, objectSegment) / ...
            (referenceNorm * objectNorm);
    end
end

[maximumScore, maximumIndex] = max(scores, [], 'omitnan');
if isempty(maximumIndex) || ~isfinite(maximumScore)
    delaySteps = NaN;
    maximumScore = NaN;
    return
end

delaySteps = maximumIndex - 1;
if maximumIndex > 1 && maximumIndex < numel(scores)
    leftScore = scores(maximumIndex-1);
    centerScore = scores(maximumIndex);
    rightScore = scores(maximumIndex+1);
    curvature = leftScore - 2*centerScore + rightScore;
    if isfinite(curvature) && curvature < 0
        fractionalOffset = 0.5 * (leftScore-rightScore) / curvature;
        fractionalOffset = max(-0.5, min(0.5, fractionalOffset));
        delaySteps = delaySteps + fractionalOffset;
    end
end
end

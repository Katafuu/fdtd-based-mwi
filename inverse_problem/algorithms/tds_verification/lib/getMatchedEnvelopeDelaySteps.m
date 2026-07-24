function [delaySteps, matchScores] = getMatchedEnvelopeDelaySteps( ...
        objectEz, referenceEz, receiverIndices, cfg, ...
        windowHalfWidthSteps, maximumDelaySteps)
%getMatchedEnvelopeDelaySteps Estimate object/reference envelope lags.

referenceArrivals = getArrivalSteps(referenceEz, receiverIndices, cfg);
numReceivers = size(receiverIndices, 1);
delaySteps = nan(numReceivers, 1);
matchScores = nan(numReceivers, 1);

for receiverIndex = 1:numReceivers
    if ~isfinite(referenceArrivals(receiverIndex))
        continue
    end

    receiverCell = receiverIndices(receiverIndex, :);
    objectSignal = reshape( ...
        objectEz(receiverCell(1), receiverCell(2), :), [], 1);
    referenceSignal = reshape( ...
        referenceEz(receiverCell(1), receiverCell(2), :), [], 1);
    [delaySteps(receiverIndex), matchScores(receiverIndex)] = ...
        estimateEnvelopeDelaySteps(referenceSignal, objectSignal, ...
        referenceArrivals(receiverIndex), windowHalfWidthSteps, ...
        maximumDelaySteps);
end
end

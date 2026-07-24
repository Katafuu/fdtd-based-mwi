function arrival_steps = getArrivalSteps(Ez, receiver_cell_indices, cfg)
%getArrivalSteps Find the first prominent received-pulse envelope peak.
%
% The active transmitter and its waveform are obtained from cfg. For each
% receiver, samples before the earliest free-space arrival are excluded.

    activeTx = find(any(cfg.source.samples ~= 0, 2));
    if numel(activeTx) ~= 1
        error('getArrivalSteps:ActiveTransmitterCount', ...
            'Exactly one row of cfg.source.samples must be active.');
    end

    sourceSignal = cfg.source.samples(activeTx, :).';
    sourceEnvelope = abs(hilbert(sourceSignal));
    sourceEnvelopeMaximum = max(sourceEnvelope);
    sourceOnsetStep = find(sourceEnvelope > 0.01*sourceEnvelopeMaximum, ...
        1, 'first');
    if isempty(sourceOnsetStep)
        error('getArrivalSteps:EmptySource', ...
            'The active source waveform must contain a nonzero pulse.');
    end

    transmitterIndex = cfg.antennas.pos(activeTx, :);
    numRx = size(receiver_cell_indices, 1);
    arrival_steps = nan(numRx, 1);

    for rx = 1:numRx
        receiverIndex = receiver_cell_indices(rx, :);
        distance = hypot( ...
            (receiverIndex(1) - transmitterIndex(1))*cfg.dx, ...
            (receiverIndex(2) - transmitterIndex(2))*cfg.dy);
        travelSteps = ceil(distance/(cfg.c0*cfg.dt));
        gateStart = max(1, sourceOnsetStep + travelSteps - 1);

        signal = reshape(Ez(receiverIndex(1), receiverIndex(2), :), [], 1);
        signalEnvelope = abs(hilbert(signal));

        preGateEnd = min(gateStart - 1, numel(signalEnvelope));
        if preGateEnd >= 1
            preGateEnvelope = signalEnvelope(1:preGateEnd);
            preGateMedian = median(preGateEnvelope);
            medianAbsoluteDeviation = median(abs(preGateEnvelope - preGateMedian));
            preGateMaximum = max(preGateEnvelope);
        else
            medianAbsoluteDeviation = 0;
            preGateMaximum = 0;
        end
        minProminence = max([1e-12, ...
            6*1.4826*medianAbsoluteDeviation, 10*preGateMaximum]);

        if gateStart > numel(signalEnvelope)
            continue;
        end

        [~, peakLocations] = findpeaks(signalEnvelope(gateStart:end), ...
            'MinPeakProminence', minProminence);
        if ~isempty(peakLocations)
            arrival_steps(rx) = gateStart + peakLocations(1) - 1;
        end
    end
end

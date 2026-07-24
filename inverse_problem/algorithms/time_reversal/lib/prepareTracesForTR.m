function [tr_data_rev, traceInfo] = prepareTracesForTR(rx_scat, dt, opts)
%prepareTracesForTR Window, normalize, and reverse scattered receiver traces.

if nargin < 3 || isempty(opts)
    opts = struct();
end
if ~isnumeric(rx_scat) || ~ismatrix(rx_scat) || any(~isfinite(rx_scat(:)))
    error('prepareTracesForTR:InvalidReceiverData', ...
        'rx_scat must be a finite numeric Nant-by-Nt array.');
end
if ~isscalar(dt) || ~isnumeric(dt) || ~isfinite(dt) || dt <= 0
    error('prepareTracesForTR:InvalidTimeStep', ...
        'dt must be a positive finite scalar.');
end

removeTxChannel = logicalScalar(getOption(opts, 'removeTxChannel', false), ...
    'prepareTracesForTR:InvalidRemoveTxChannel');
txAntenna = getOption(opts, 'txAntenna', []);
applyTemporalWindow = logicalScalar(getOption(opts, 'applyTemporalWindow', true), ...
    'prepareTracesForTR:InvalidApplyTemporalWindow');
normalizeTraces = logicalScalar(getOption(opts, 'normalizeTraces', false), ...
    'prepareTracesForTR:InvalidNormalizeTraces');
temporalWindowTau = getOption(opts, 'temporalWindowTau', 25e-12);
if ~isscalar(temporalWindowTau) || ~isnumeric(temporalWindowTau) || ...
        ~isfinite(temporalWindowTau) || temporalWindowTau <= 0
    error('prepareTracesForTR:InvalidTemporalWindowTau', ...
        'temporalWindowTau must be a positive finite scalar.');
end

[numReceivers, numSamples] = size(rx_scat);
trData = rx_scat;
if removeTxChannel
    if isempty(txAntenna) || ~isscalar(txAntenna) || txAntenna ~= round(txAntenna) || ...
            txAntenna < 1 || txAntenna > numReceivers
        error('prepareTracesForTR:InvalidTxAntenna', ...
            'opts.txAntenna must index one receiver when removeTxChannel is true.');
    end
    trData(txAntenna, :) = 0;
end

trDataUnwindowed = trData;
time = (0:numSamples-1) .* dt;
temporalWindow = ones(numReceivers, numSamples);
peakIndex = ones(numReceivers, 1);
peakTime = zeros(numReceivers, 1);

if applyTemporalWindow
    for rxIdx = 1:numReceivers
        if any(trData(rxIdx, :) ~= 0)
            [~, peakIndex(rxIdx)] = max(abs(trData(rxIdx, :)));
            peakTime(rxIdx) = time(peakIndex(rxIdx));
            temporalWindow(rxIdx, :) = ...
                exp(-((time - peakTime(rxIdx)) ./ temporalWindowTau).^2);
            trData(rxIdx, :) = trData(rxIdx, :) .* temporalWindow(rxIdx, :);
        end
    end
end

trDataWindowed = trData;
tr_data_rev = fliplr(trDataWindowed);
traceScale = ones(numReceivers, 1);
if normalizeTraces
    traceScale = max(abs(tr_data_rev), [], 2);
    traceScale(traceScale == 0) = 1;
    tr_data_rev = tr_data_rev ./ traceScale;
end

traceInfo = struct();
traceInfo.raw = rx_scat;
traceInfo.unwindowed = trDataUnwindowed;
traceInfo.windowed = trDataWindowed;
traceInfo.reversed = tr_data_rev;
traceInfo.temporalWindow = temporalWindow;
traceInfo.peakIndex = peakIndex;
traceInfo.peakTime = peakTime;
traceInfo.time = time;
traceInfo.traceScale = traceScale;
traceInfo.removeTxChannel = removeTxChannel;
traceInfo.applyTemporalWindow = applyTemporalWindow;
traceInfo.temporalWindowTau = temporalWindowTau;
traceInfo.normalizeTraces = normalizeTraces;
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end

function value = logicalScalar(value, errorId)
if ~isscalar(value) || ~(islogical(value) || isnumeric(value))
    error(errorId, 'Expected a logical scalar.');
end
value = logical(value);
end

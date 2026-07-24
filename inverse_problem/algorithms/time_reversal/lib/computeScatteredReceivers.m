function rx_scat = computeScatteredReceivers(result_tot, result_inc)
%computeScatteredReceivers Subtract incident receiver traces from total traces.

if ~isstruct(result_tot) || ~isfield(result_tot, 'rx_signals')
    error('computeScatteredReceivers:MissingTotalReceivers', ...
        'result_tot.rx_signals is required.');
end
if ~isstruct(result_inc) || ~isfield(result_inc, 'rx_signals')
    error('computeScatteredReceivers:MissingIncidentReceivers', ...
        'result_inc.rx_signals is required.');
end
if ~isequal(size(result_tot.rx_signals), size(result_inc.rx_signals))
    error('computeScatteredReceivers:ReceiverSizeMismatch', ...
        'Total and incident receiver arrays must have matching sizes.');
end

rx_scat = result_tot.rx_signals - result_inc.rx_signals;
end

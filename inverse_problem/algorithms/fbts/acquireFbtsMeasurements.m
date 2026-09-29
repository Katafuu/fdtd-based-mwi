function measurements = acquireFbtsMeasurements(cfg)
%acquireFbtsMeasurements Acquire configured TXs and every configured receiver.
[cfg, time, pulse] = prepareFbtsSource(cfg);
tx = cfg.antennas.txAntennas(:).';
n = cfg.antennas.numAntennas;
validateattributes(tx, {'numeric'}, {'vector','integer','>=',1,'<=',n});
assert(numel(unique(tx)) == numel(tx), 'fbts:InvalidTransmitters', 'Duplicate transmitters.');
original = 1:n;
if isfield(cfg.antennas, 'originalIndices'), original = cfg.antennas.originalIndices; end
rx = zeros(numel(tx), n, cfg.Nt);
seconds = zeros(numel(tx), 1);
startedUtc=fbtsUtcNow(); started = tic;
for j = 1:numel(tx)
    c = cfg;
    c.returnEz = false; c.returnHx = false; c.returnHy = false;
    c.returnRxSignals = true;
    c.source.samples = zeros(n, cfg.Nt);
    c.source.samples(tx(j), :) = pulse;
    timer = tic;
    out = fbtsSolve(c);
    seconds(j) = toc(timer);
    rx(j,:,:) = reshape(out.rx_signals, 1, n, cfg.Nt);
end
identity = fbtsAcquisitionIdentity(cfg);
assert(all(isfinite(rx),'all'),'fbts:NonfiniteMeasurements','The FDTD acquisition contains nonfinite samples.');
metadata = struct('schema_version', 1, 'identity', identity, ...
    'tx_ids', original(tx), 'rx_ids', original, ...
    'tx_positions', cfg.antennas.pos(tx,:), 'rx_positions', cfg.antennas.pos, ...
    'time', time, 'dimension_order', 'tx,rx,time', 'complete', true, ...
    'noise', 'none', 'signal', 'total Ez', 'self_receiver_included', true);
metadata.rx_hash = fbtsHash(rx);
metadata.fingerprint = fbtsHash(metadata);
measurements = struct('rx', rx, 'metadata', metadata, ...
    'timing', struct('started_utc',startedUtc,'finished_utc',fbtsUtcNow(), ...
    'solver_seconds', seconds, 'wall_seconds', toc(started), ...
    'solver_count', numel(tx)));
end

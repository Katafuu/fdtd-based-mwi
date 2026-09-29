function cfg = subsetFbtsConfig(baseCfg, active)
%subsetFbtsConfig Active antennas in ascending original order; all TX and RX.
active = sort(active(:).');
n = baseCfg.antennas.numAntennas;
validateattributes(active, {'numeric'}, {'nonempty','integer','>=',1,'<=',n});
assert(numel(unique(active)) == numel(active), 'fbts:InvalidSubset', 'Duplicate active indices.');
cfg = baseCfg;
cfg.antennas.pos = baseCfg.antennas.pos(active,:);
cfg.antennas.numAntennas = numel(active);
cfg.antennas.txAntennas = 1:numel(active);
cfg.antennas.originalIndices = active;
cfg.source.location = cfg.antennas.pos;
cfg.source.samples = zeros(numel(active), cfg.Nt);
end

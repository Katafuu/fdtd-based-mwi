function data = loadFbtsBatchCase(resultFile, verify)
%loadFbtsBatchCase Assemble a portable case and the runnable reduced cfg.
if nargin < 2, verify = true; end
result = load(resultFile);
root = fullfile(fileparts(resultFile),result.caseInfo.dataset_root_relative);
refs = result.caseInfo.references;
paths = {fullfile(root,refs.system_file),fullfile(root,refs.truth_file), ...
    fullfile(root,refs.acquisition_file)};
if verify
    hashes = {refs.system_hash,refs.truth_hash,refs.acquisition_hash};
    for k = 1:3
        assert(strcmp(fbtsHash(paths{k},'file'),hashes{k}), ...
            'fbts:DatasetMismatch','Shared file missing or modified: %s.',paths{k});
    end
    status = load(fullfile(fileparts(resultFile),'status.mat'),'status');
    assert(status.status.status == "success" && strcmp(status.status.result_hash,fbtsHash(resultFile,'file')), ...
        'fbts:DatasetMismatch','Result is incomplete or modified.');
end
system = load(paths{1}); truth = load(paths{2}); acquisition = load(paths{3});
cfg = system.model.cfg_template;
g = system.model.grid_metadata;
cfg.grid = g;
[cfg.grid.xIndex,cfg.grid.yIndex] = ndgrid(1:cfg.Nx,1:cfg.Ny);
cfg.grid.xPhysical = g.originPhysical(1)+(cfg.grid.xIndex-1)*cfg.dx;
cfg.grid.yPhysical = g.originPhysical(2)+(cfg.grid.yIndex-1)*cfg.dy;
for key = {'epsr','cond_e','cond_m','murx','mury'}
    cfg.grid.(key{1}) = truth.([key{1} '_true']);
    cfg.grid.([key{1} '_bg']) = system.([key{1} '_background']);
end
cfg.targets = truth.targets;
for k = 1:numel(cfg.targets), cfg.targets(k).mask = truth.target_masks(:,:,k); end
cfg.antennas.doiMask = system.doi_mask;
cfg.source.samples = zeros(cfg.antennas.numAntennas,cfg.Nt);
cfg = subsetFbtsConfig(cfg,result.caseInfo.entry.active_indices);
[cfg,~,~,~] = prepareFbtsSource(cfg);
measurements = struct('rx',acquisition.rx_full,'metadata',acquisition.measurement, ...
    'timing',acquisition.measurementTiming);
% Analysis must not require a MEX binary, let alone the original platform's
% binary. runFbts validates physics/solver identity when actually rerunning.
active = result.caseInfo.tx_order_original;
[hasTx,tx] = ismember(active,acquisition.measurement.tx_ids);
[hasRx,rx] = ismember(active,acquisition.measurement.rx_ids);
assert(all(hasTx) && all(hasRx),'fbts:DatasetMismatch','Case antenna IDs are missing from the saved scan.');
data = struct('result',result,'system',system,'truth',truth, ...
    'acquisition',acquisition,'cfg',cfg,'measurements',measurements, ...
    'rx_selected',acquisition.rx_full(tx,rx,:));
end

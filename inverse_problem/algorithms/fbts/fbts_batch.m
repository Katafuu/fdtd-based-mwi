% fbts_batch Run reflection/rotation-reduced FBTS batches with shared RX scans.
% Build cfg first, or provide fbtsScenes (cell/struct array of full scenes).
% Optional: fbtsOptions (algorithm/output), fbtsBatchOptions (batch controls).
if exist('fbtsScenes','var')
    activeFbtsScenes = fbtsScenes;
else
    assert(exist('cfg','var') == 1 && isstruct(cfg), ...
        'fbts:MissingConfig','Run build_cfg.m or provide fbtsScenes first.');
    activeFbtsScenes = cfg;
end
if exist('fbtsOptions','var'), activeFbtsOptions = resolveFbtsOptions(fbtsOptions);
else, activeFbtsOptions = resolveFbtsOptions(); end
if exist('fbtsBatchOptions','var'), activeFbtsBatchOptions = resolveFbtsBatchOptions(fbtsBatchOptions);
else, activeFbtsBatchOptions = resolveFbtsBatchOptions(); end
[batchSummary,batchOutputDirectory,batchManifest] = runFbtsBatch( ...
    activeFbtsScenes,activeFbtsOptions,activeFbtsBatchOptions);

function coarseCfg = prepareCoarseSensitivityCfg(cfg, downsampleFactor)
%prepareCoarseSensitivityCfg Build an executable coarse FBTS sensitivity cfg.
%   Spatial resampling deliberately leaves temporal fields unchanged, so this
%   helper immediately recomputes the CFL timestep, the minimum Nt implied by
%   cfg.deltaF, and the source sample matrix required by fdtd_mex.

if ~isfield(cfg, 'deltaF') || ~isnumeric(cfg.deltaF) || ...
        ~isscalar(cfg.deltaF) || ~isfinite(cfg.deltaF) || cfg.deltaF <= 0
    error('fbts:prepareCoarseSensitivityCfg:InvalidDeltaF', ...
        'cfg.deltaF must be a finite positive scalar.');
end

coarseCfg = fdtdmat.downsampleCfg(cfg, downsampleFactor);
coarseCfg.dt = 1 / (coarseCfg.c0 * ...
    sqrt(1/coarseCfg.dx^2 + 1/coarseCfg.dy^2));
coarseCfg.Nt = ceil(1 / (coarseCfg.dt * coarseCfg.deltaF));
coarseCfg.snapshotStart = 0;
coarseCfg.snapshotStride = 0;
coarseCfg.returnEz = false;
coarseCfg.returnHx = false;
coarseCfg.returnHy = false;
coarseCfg.returnRxSignals = true;
coarseCfg.solverThreads = 1;
coarseCfg.solverBackend = 'cpu';
coarseCfg.source.time = (0:coarseCfg.Nt-1) .* coarseCfg.dt;
coarseCfg.source.location = coarseCfg.antennas.pos;
coarseCfg.source.samples = zeros( ...
    coarseCfg.antennas.numAntennas, coarseCfg.Nt);
end

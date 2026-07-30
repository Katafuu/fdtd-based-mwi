function resizedCfg = upsampleCfg(cfg, factor)
%upsampleCfg Spatially upsample an FDTD configuration by an integer factor.
%   resizedCfg = fdtdmat.upsampleCfg(cfg, factor) preserves the physical x/y
%   size while increasing Nx and Ny and reducing dx and dy. Grid material
%   maps, PML maps, antennas, DOI/target masks, and index geometry are updated.
%
%   WARNING: This is deliberately a spatial-only transform. It does not
%   change cfg.dt, cfg.Nt, cfg.source.time, or cfg.source.samples. The old dt
%   will generally violate the CFL limit on the finer grid. Before executing
%   the returned configuration, recompute dt and Nt and rebuild time samples.

resizedCfg = fdtdmat.utility_resampleCfgSpatial(cfg, factor, "up");
end

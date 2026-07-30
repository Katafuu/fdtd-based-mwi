function resizedCfg = downsampleCfg(cfg, factor)
%downsampleCfg Spatially downsample an FDTD configuration by an integer factor.
%   resizedCfg = fdtdmat.downsampleCfg(cfg, factor) preserves the physical
%   x/y size while reducing Nx and Ny and increasing dx and dy. Grid material
%   maps, PML maps, antennas, DOI/target masks, and index geometry are updated.
%
%   WARNING: This is deliberately a spatial-only transform. It does not
%   change cfg.dt, cfg.Nt, cfg.source.time, or cfg.source.samples. Before the
%   returned configuration is executed, recompute a stable dt, recompute Nt
%   (for example, ceil(1/(dt*deltaF))), and rebuild all time-sampled sources.

resizedCfg = fdtdmat.utility_resampleCfgSpatial(cfg, factor, "down");
end

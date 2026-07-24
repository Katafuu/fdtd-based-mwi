function cfg = setBackgroundDefault(cfg, background)
%setBackgroundDefault Replace cfg.grid with a canonical background grid.

originPhysical = [0 0];
if isfield(cfg, 'grid') && isstruct(cfg.grid) && ...
        isfield(cfg.grid, 'originPhysical') && ~isempty(cfg.grid.originPhysical)
    originPhysical = cfg.grid.originPhysical;
end
cfg.grid = fdtdmat.createGrid([cfg.Nx cfg.Ny], [cfg.dx cfg.dy], ...
    background, originPhysical);
end

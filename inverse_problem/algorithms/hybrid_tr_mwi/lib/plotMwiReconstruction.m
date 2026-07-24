function imageHandle = plotMwiReconstruction(axesHandle, cfg, epsrAvgFull)
%plotMwiReconstruction Plot a full-grid MWI permittivity reconstruction.

if ~isequal(size(epsrAvgFull), [cfg.Nx cfg.Ny])
    error('plotMwiReconstruction:GridSizeMismatch', ...
        'epsrAvgFull must have size [cfg.Nx cfg.Ny].');
end

imageHandle = imagesc(axesHandle, 1:cfg.Nx, 1:cfg.Ny, epsrAvgFull.', ...
    'Tag', 'epsrReconstruction');
axis(axesHandle, 'equal');
xlim(axesHandle, [0.5 cfg.Nx + 0.5]);
ylim(axesHandle, [0.5 cfg.Ny + 0.5]);
set(axesHandle, 'YDir', 'normal');
xlabel(axesHandle, 'x index');
ylabel(axesHandle, 'y index');
colorbar(axesHandle);
title(axesHandle, ...
    'Circular-array weighted average \epsilon_r reconstruction, full grid');

plotDomainOverlays(axesHandle, cfg, cfg.antennas.txAntennas(1), 'w');
end

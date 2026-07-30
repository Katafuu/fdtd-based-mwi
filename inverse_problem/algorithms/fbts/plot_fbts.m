% plot_fbts Display the true and reconstructed relative permittivity.

colorLimits = [ ...
    min(results.epsr_true(:)), ...
    max(results.epsr_true(:))];

trueFigure = figure('Name', 'True relative permittivity', ...
    'NumberTitle', 'off');
trueAxes = axes('Parent', trueFigure);
imagesc(trueAxes, 1:cfg.Nx, 1:cfg.Ny, results.epsr_true.');
axis(trueAxes, 'equal');
xlim(trueAxes, [0.5 cfg.Nx + 0.5]);
ylim(trueAxes, [0.5 cfg.Ny + 0.5]);
set(trueAxes, 'YDir', 'normal');
xlabel(trueAxes, 'x index');
ylabel(trueAxes, 'y index');
clim(trueAxes, colorLimits);
colorbar(trueAxes);
title(trueAxes, 'True relative permittivity, full grid');
plotDomainOverlays(trueAxes, cfg, [], 'w');

estimatedFigure = figure('Name', 'Estimated relative permittivity', ...
    'NumberTitle', 'off');
estimatedAxes = axes('Parent', estimatedFigure);
imagesc(estimatedAxes, 1:cfg.Nx, 1:cfg.Ny, results.epsr_est.');
axis(estimatedAxes, 'equal');
xlim(estimatedAxes, [0.5 cfg.Nx + 0.5]);
ylim(estimatedAxes, [0.5 cfg.Ny + 0.5]);
set(estimatedAxes, 'YDir', 'normal');
xlabel(estimatedAxes, 'x index');
ylabel(estimatedAxes, 'y index');
clim(estimatedAxes, colorLimits);
colorbar(estimatedAxes);
title(estimatedAxes, 'FBTS relative permittivity reconstruction, %d iterations', numIterations);
plotDomainOverlays(estimatedAxes, cfg, [], 'w');

convergenceFigure = figure('Name', 'FBTS convergence', ...
    'NumberTitle', 'off');
convergenceAxes = axes('Parent', convergenceFigure);
plot(convergenceAxes, results.iteration, ...
    results.relative_error_percent, 'o-', 'LineWidth', 1.2);
grid(convergenceAxes, 'on');
xlabel(convergenceAxes, 'FBTS iteration');
ylabel(convergenceAxes, 'Relative permittivity error (%)');
title(convergenceAxes, 'FBTS DOI relative-error convergence');

costFigure = figure('Name', 'FBTS cost convergence', ...
    'NumberTitle', 'off');
costAxes = axes('Parent', costFigure);
plot(costAxes, results.iteration, results.cost, ...
    'o-', 'LineWidth', 1.2);
grid(costAxes, 'on');
xlabel(costAxes, 'FBTS iteration');
ylabel(costAxes, 'Weighted data cost');
title(costAxes, 'FBTS cost convergence');

figureFiles = results.output_files.figures;
exportgraphics(trueFigure, figureFiles(1), 'Resolution', 300);
exportgraphics(estimatedFigure, figureFiles(2), 'Resolution', 300);
exportgraphics(convergenceFigure, figureFiles(3), 'Resolution', 300);
exportgraphics(costFigure, figureFiles(4), 'Resolution', 300);

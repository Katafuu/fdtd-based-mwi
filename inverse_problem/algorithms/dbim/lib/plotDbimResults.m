function figureFiles = plotDbimResults( ...
    trueCfg, dbimResult, outputDirectory, showPlots, resolution)
%plotDbimResults Plot and save DBIM maps and convergence histories.

if ~isstruct(trueCfg) || ~isfield(trueCfg, 'grid') || ...
        ~isstruct(dbimResult) || ~isfield(dbimResult, 'estimate') || ...
        ~isfield(dbimResult, 'history')
    error('plotDbimResults:InvalidInput', ...
        'Expected a true cfg and a DBIM result structure.');
end
if ~(ischar(outputDirectory) || ...
        (isstring(outputDirectory) && isscalar(outputDirectory))) || ...
        strlength(string(outputDirectory)) == 0
    error('plotDbimResults:InvalidOutputDirectory', ...
        'outputDirectory must be a nonempty text scalar.');
end
if ~isscalar(showPlots) || ...
        ~(islogical(showPlots) || isnumeric(showPlots))
    error('plotDbimResults:InvalidShowPlots', ...
        'showPlots must be a logical scalar.');
end
if ~isnumeric(resolution) || ~isscalar(resolution) || ...
        ~isfinite(resolution) || resolution <= 0
    error('plotDbimResults:InvalidResolution', ...
        'resolution must be a finite positive scalar.');
end

outputDirectory = char(outputDirectory);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
figureVisibility = 'off';
if logical(showPlots)
    figureVisibility = 'on';
end

x = (0:trueCfg.Nx-1) .* trueCfg.dx;
y = (0:trueCfg.Ny-1) .* trueCfg.dy;

materialFigure = figure( ...
    'Name', 'DBIM material reconstruction', ...
    'Visible', figureVisibility);
layout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plotMaterial(x, y, trueCfg.grid.epsr, 'True relative permittivity', '\epsilon_r');
nexttile;
plotMaterial(x, y, dbimResult.estimate.epsr, ...
    'Reconstructed relative permittivity', '\epsilon_r');
nexttile;
plotMaterial(x, y, trueCfg.grid.cond_e, ...
    'True conductivity', '\sigma [S/m]');
nexttile;
plotMaterial(x, y, dbimResult.estimate.sigma, ...
    'Reconstructed conductivity', '\sigma [S/m]');
title(layout, sprintf('FDTD-DBIM result (best iteration %d)', ...
    dbimResult.bestIteration));

convergenceFigure = figure( ...
    'Name', 'DBIM convergence', ...
    'Visible', figureVisibility);
semilogy(dbimResult.history.evaluatedIteration, ...
    dbimResult.history.residual, 'o-', 'LineWidth', 1.2, ...
    'DisplayName', 'Normalized data residual');
hold on;
validUpdate = isfinite(dbimResult.history.updateNorm);
semilogy(dbimResult.history.evaluatedIteration(validUpdate), ...
    dbimResult.history.updateNorm(validUpdate), 's-', 'LineWidth', 1.2, ...
    'DisplayName', 'Relative update norm');
hold off;
grid on;
xlabel('Applied DBIM updates');
ylabel('Relative norm');
title('DBIM convergence');
legend('Location', 'best');

figureFiles = [
    string(fullfile(outputDirectory, '01_dbim_material_reconstruction.png'))
    string(fullfile(outputDirectory, '02_dbim_convergence.png'))
];
exportgraphics(materialFigure, figureFiles(1), 'Resolution', resolution);
exportgraphics(convergenceFigure, figureFiles(2), 'Resolution', resolution);

if ~logical(showPlots)
    close([materialFigure convergenceFigure]);
end
end

function plotMaterial(x, y, materialMap, plotTitle, colorbarLabel)
imagesc(x, y, double(materialMap).');
axis xy image;
xlabel('x [m]');
ylabel('y [m]');
title(plotTitle);
colorbarHandle = colorbar;
colorbarHandle.Label.String = colorbarLabel;
end
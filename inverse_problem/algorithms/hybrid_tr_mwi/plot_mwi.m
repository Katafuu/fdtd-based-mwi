function plot_mwi(cfg, plotData, outputDir)
%plot_mwi Create and export the six hybrid TR/MWI diagnostic figures.

requiredFields = {'trResult', 'focusPoint', 'mask4D', 'xRange', ...
    'yRange', 'pairAccepted', 'epsrAvgFull'};
missingFields = requiredFields(~isfield(plotData, requiredFields));
if ~isempty(missingFields)
    error('plot_mwi:MissingPlotData', ...
        'plotData is missing required field(s): %s.', ...
        strjoin(missingFields, ', '));
end

if ~isfolder(outputDir)
    mkdir(outputDir);
end

figuresBefore = findall(groot, 'Type', 'figure');
previousDefaultVisibility = get(groot, 'DefaultFigureVisible');
plotCleanup = onCleanup(@() restorePlotState( ...
    figuresBefore, previousDefaultVisibility));
set(groot, 'DefaultFigureVisible', 'off');

antennaFigure = figure('Visible', 'off');
antennaAxes = axes('Parent', antennaFigure);
axes(antennaAxes);
plotAntennaArray(cfg.antennas.pos, cfg.antennas.center, ...
    cfg.Nx, cfg.Ny, cfg.pml.thickness);
title(antennaAxes, 'Circular antenna array');

trFigures = plotTrComparison(cfg, plotData.trResult);
maxFocusPosition = get(trFigures(2), 'Position');

footprintFigure = figure( ...
    'Name', 'Accepted TR-shifted footprint coverage', ...
    'NumberTitle', 'off', 'Visible', 'off', ...
    'Position', maxFocusPosition);
footprintAxes = axes('Parent', footprintFigure);
axes(footprintAxes);
plotCircularFootprintCoverage(plotData.mask4D, plotData.xRange, ...
    plotData.yRange, plotData.pairAccepted, cfg, plotData.focusPoint);

reconstructionFigure = figure('Visible', 'off', ...
    'Position', maxFocusPosition);
reconstructionAxes = axes('Parent', reconstructionFigure);
plotMwiReconstruction(reconstructionAxes, cfg, plotData.epsrAvgFull);
drawnow;

maxFocusFile = fullfile(outputDir, '03_tr_max_magnitude_focus.png');
exportgraphics(trFigures(2), maxFocusFile, 'Resolution', 300);
maxFocusInfo = imfinfo(maxFocusFile);

exportgraphics(antennaFigure, ...
    fullfile(outputDir, '01_antenna_array.png'), 'Resolution', 300);
exportgraphics(trFigures(1), ...
    fullfile(outputDir, '02_tr_actual_target.png'), 'Resolution', 300);
exportgraphics(trFigures(3), ...
    fullfile(outputDir, '04_tr_minimum_r_focus.png'), 'Resolution', 300);

exportMatchedPng(footprintAxes, ...
    fullfile(outputDir, '05_footprint_coverage.png'), ...
    maxFocusInfo.Width, maxFocusInfo.Height);
exportMatchedPng(reconstructionAxes, ...
    fullfile(outputDir, '06_epsr_reconstruction.png'), ...
    maxFocusInfo.Width, maxFocusInfo.Height);
end

function exportMatchedPng(graphicsHandle, filename, targetWidth, targetHeight)
exportgraphics(graphicsHandle, filename, 'Resolution', 300, ...
    'Width', targetWidth, 'Height', targetHeight, ...
    'Units', 'pixels', 'Padding', 0, 'PreserveAspectRatio', 'off');
exportedInfo = imfinfo(filename);
if exportedInfo.Width ~= targetWidth || exportedInfo.Height ~= targetHeight
    error('plot_mwi:ExportSizeMismatch', ...
        'Could not export %s at exactly %d-by-%d pixels.', ...
        filename, targetWidth, targetHeight);
end
end

function restorePlotState(figuresBefore, previousDefaultVisibility)
set(groot, 'DefaultFigureVisible', previousDefaultVisibility);
createdFigures = setdiff(findall(groot, 'Type', 'figure'), ...
    figuresBefore, 'stable');
close(createdFigures(isgraphics(createdFigures, 'figure')));
end

function result = tds_2d(cfg)
%tds_2d Reconstruct one effective epsr inside time-reversal support.
% With no input, build the default hybrid configuration. Both call forms
% save shape figures, truth, reconstruction, and DOI error in a new run folder.

if nargin < 1 || isempty(cfg)
    cfg = build_cfg(struct());
end
assert(isstruct(cfg) && isscalar(cfg), ...
    'tds_2d:InvalidConfig', 'cfg must be a scalar struct.');

outputRoot = fullfile(fileparts(mfilename('fullpath')), 'figs', 'tds_2d');
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end
existingRuns = dir(fullfile(outputRoot, 'run_*'));
runNumbers = zeros(0, 1);
for runIdx = 1:numel(existingRuns)
    if ~existingRuns(runIdx).isdir
        continue;
    end
    token = regexp(existingRuns(runIdx).name, ...
        '^run_(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        runNumbers(end + 1, 1) = str2double(token{1}); %#ok<AGROW>
    end
end
outputDir = fullfile(outputRoot, ...
    sprintf('run_%04d', max([0; runNumbers]) + 1));
mkdir(outputDir);

[supportMask, shapeDetails] = shape_estimate( ...
    cfg, fullfile(outputDir, 'shape_estimate'));
measurement = acquireTransmissionMeasurements(cfg);
[recoveredTargetEpsr, maximumChordDetails] = ...
    estimateMaximumChordEpsr(measurement, supportMask);
reconstructedEpsr = reconstructHomogeneousTargetMap( ...
    supportMask, measurement.backgroundEpsr, recoveredTargetEpsr);

trueEpsr = cfg.grid.epsr;
doiMask = logical(cfg.antennas.doiMask);
assert(isequal(size(trueEpsr), size(reconstructedEpsr), size(doiMask)) ...
    && any(doiMask(:)), 'tds_2d:InvalidDoi', ...
    'Truth, reconstruction, and nonempty DOI mask must share a grid size.');
doiDifference = reconstructedEpsr(doiMask) - trueEpsr(doiMask);
doiError = struct( ...
    'numPixels', nnz(doiMask), ...
    'meanAbsoluteError', mean(abs(doiDifference)), ...
    'rootMeanSquareError', sqrt(mean(doiDifference.^2)), ...
    'meanAbsolutePercentageError', ...
        100 * mean(abs(doiDifference) ./ abs(trueEpsr(doiMask))));

result = struct( ...
    'outputDir', outputDir, ...
    'supportMask', supportMask, ...
    'recoveredTargetEpsr', recoveredTargetEpsr, ...
    'reconstructedEpsr', reconstructedEpsr, ...
    'trueEpsr', trueEpsr, ...
    'doiError', doiError, ...
    'shapeDetails', shapeDetails, ...
    'maximumChordDetails', maximumChordDetails, ...
    'measurement', measurement);

reconstructionFigure = figure('Visible', 'off');
figureCleanup = onCleanup(@() close(reconstructionFigure)); %#ok<NASGU>
axesHandle = axes('Parent', reconstructionFigure);
imagesc(axesHandle, 1:cfg.Nx, 1:cfg.Ny, reconstructedEpsr.');
axis(axesHandle, 'equal', 'tight');
set(axesHandle, 'YDir', 'normal');
colormap(axesHandle, turbo);
colorLimits = [min([trueEpsr(:); reconstructedEpsr(:)]), ...
    max([trueEpsr(:); reconstructedEpsr(:)])];
if colorLimits(1) < colorLimits(2)
    clim(axesHandle, colorLimits);
end
colorbar(axesHandle);
xlabel(axesHandle, 'x grid index');
ylabel(axesHandle, 'y grid index');
title(axesHandle, sprintf( ...
    'Estimated support: reconstructed \\epsilon_r = %.4f', ...
    recoveredTargetEpsr));
exportgraphics(reconstructionFigure, ...
    fullfile(outputDir, 'epsr_reconstruction.png'), 'Resolution', 300);

truthFigure = figure('Visible', 'off');
truthCleanup = onCleanup(@() close(truthFigure)); %#ok<NASGU>
truthAxes = axes('Parent', truthFigure);
imagesc(truthAxes, 1:cfg.Nx, 1:cfg.Ny, trueEpsr.');
axis(truthAxes, 'equal', 'tight');
set(truthAxes, 'YDir', 'normal');
colormap(truthAxes, turbo);
if colorLimits(1) < colorLimits(2)
    clim(truthAxes, colorLimits);
end
colorbar(truthAxes);
xlabel(truthAxes, 'x grid index');
ylabel(truthAxes, 'y grid index');
title(truthAxes, 'True relative permittivity');
exportgraphics(truthFigure, ...
    fullfile(outputDir, 'epsr_true.png'), 'Resolution', 300);

errorFile = fopen(fullfile(outputDir, 'doi_error.txt'), 'w');
assert(errorFile ~= -1, 'tds_2d:ErrorFileOpen', ...
    'Could not create DOI error file in %s.', outputDir);
errorFileCleanup = onCleanup(@() fclose(errorFile)); %#ok<NASGU>
fprintf(errorFile, 'DOI = cfg.antennas.doiMask (%d grid cells)\n', ...
    doiError.numPixels);
fprintf(errorFile, 'MAE = mean(abs(reconstructed epsr - true epsr)) = %.9g\n', ...
    doiError.meanAbsoluteError);
fprintf(errorFile, 'RMSE = sqrt(mean((reconstructed epsr - true epsr)^2)) = %.9g\n', ...
    doiError.rootMeanSquareError);
fprintf(errorFile, 'MAPE = mean(abs(reconstructed epsr - true epsr) / abs(true epsr)) * 100 = %.9g %%\n', ...
    doiError.meanAbsolutePercentageError);
save(fullfile(outputDir, 'reconstruction.mat'), 'result');

fprintf('Estimated target epsr: %.6f from %d directed chords.\n', ...
    recoveredTargetEpsr, maximumChordDetails.numSelectedDirectedRays);
fprintf('Saved reconstruction in %s.\n', outputDir);
end

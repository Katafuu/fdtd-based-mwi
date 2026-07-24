% tds_2d_improved Reconstruct a homogeneous target with the validated method.
% This script is standalone. Set verificationSeed and forceRecompute first
% to override their defaults.

if ~exist('verificationSeed', 'var') || isempty(verificationSeed)
    verificationSeed = 1;
end
if ~exist('forceRecompute', 'var') || isempty(forceRecompute)
    forceRecompute = false;
end

measurement = runMeasurementCase(verificationSeed, logical(forceRecompute));
verificationResults = analyzeVerificationCases({measurement});
[recoveredTargetEpsr, maximumChordDetails] = ...
    estimateMaximumChordEpsr(measurement);

backgroundEpsr = measurement.backgroundEpsr;
knownTargetEpsr = measurement.targetEpsr;
targetMask = logical(measurement.targetMask);
baselineTargetEpsr = verificationResults.baselineEpsr(1);
baselineErrorPercent = verificationResults.baselineErrorPercent(1);
improvedErrorPercent = verificationResults.maximumChordErrorPercent(1);

groundTruth = reconstructHomogeneousTargetMap( ...
    targetMask, backgroundEpsr, knownTargetEpsr);
baselineReconstruction = reconstructHomogeneousTargetMap( ...
    targetMask, backgroundEpsr, baselineTargetEpsr);
improvedReconstruction = reconstructHomogeneousTargetMap( ...
    targetMask, backgroundEpsr, recoveredTargetEpsr);

fprintf('\nValidated maximum-chord transmission-delay reconstruction\n');
fprintf('  Seed:                         %d\n', verificationSeed);
fprintf('  Known target epsr:            %.6f\n', knownTargetEpsr);
fprintf('  Original 15-pixel estimate:   %.6f (error %.2f%%)\n', ...
    baselineTargetEpsr, baselineErrorPercent);
fprintf('  Maximum-chord estimate:       %.6f (error %.2f%%)\n', ...
    recoveredTargetEpsr, improvedErrorPercent);
fprintf('  Maximum target chord:         %.3f mm\n', ...
    1e3*maximumChordDetails.maximumChordLength);
fprintf('  Selected directed rays:       %d\n\n', ...
    maximumChordDetails.numSelectedDirectedRays);

%% Known mask and selected maximum-chord paths
rayFigure = figure('Name', 'Selected maximum-chord rays', ...
    'NumberTitle', 'off');
rayAxes = axes('Parent', rayFigure);
imagesc(rayAxes, 1:measurement.gridSize(1), ...
    1:measurement.gridSize(2), double(targetMask.'));
axis(rayAxes, 'equal');
xlim(rayAxes, [0.5 measurement.gridSize(1)+0.5]);
ylim(rayAxes, [0.5 measurement.gridSize(2)+0.5]);
set(rayAxes, 'YDir', 'normal');
colormap(rayAxes, gray);
xlabel(rayAxes, 'x index');
ylabel(rayAxes, 'y index');
title(rayAxes, sprintf( ...
    'Known mask and maximum-chord paths (%.3f mm)', ...
    1e3*maximumChordDetails.maximumChordLength));
hold(rayAxes, 'on');
selectedPairs = maximumChordDetails.selectedPairs;
numAntennas = size(measurement.antennaPositions, 1);
for tx = 1:numAntennas-1
    for rx = tx+1:numAntennas
        if selectedPairs(tx, rx) || selectedPairs(rx, tx)
            endpoints = measurement.antennaPositions([tx rx], :);
            plot(rayAxes, endpoints(:, 1), endpoints(:, 2), ...
                'r-', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
end
plot(rayAxes, measurement.antennaPositions(:, 1), ...
    measurement.antennaPositions(:, 2), 'ko', ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
    'DisplayName', 'Antennas');
hold(rayAxes, 'off');

%% Ground truth, original estimate, and improved estimate
comparisonFigure = figure('Name', 'Maximum-chord quantitative verification', ...
    'NumberTitle', 'off');
comparisonLayout = tiledlayout(comparisonFigure, 1, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
colorLimits = [min([backgroundEpsr knownTargetEpsr ...
    baselineTargetEpsr recoveredTargetEpsr]), ...
    max([backgroundEpsr knownTargetEpsr ...
    baselineTargetEpsr recoveredTargetEpsr])];
if colorLimits(1) == colorLimits(2)
    colorLimits = colorLimits + [-0.5 0.5];
end

truthAxes = nexttile(comparisonLayout);
imagesc(truthAxes, groundTruth.');
formatMapAxes(truthAxes, measurement.gridSize, colorLimits, ...
    sprintf('Truth: \\epsilon_r = %.4f', knownTargetEpsr));

baselineAxes = nexttile(comparisonLayout);
imagesc(baselineAxes, baselineReconstruction.');
formatMapAxes(baselineAxes, measurement.gridSize, colorLimits, ...
    sprintf('Original: %.4f (%.2f%% error)', ...
    baselineTargetEpsr, baselineErrorPercent));

improvedAxes = nexttile(comparisonLayout);
imagesc(improvedAxes, improvedReconstruction.');
formatMapAxes(improvedAxes, measurement.gridSize, colorLimits, ...
    sprintf('Maximum chord: %.4f (%.2f%% error)', ...
    recoveredTargetEpsr, improvedErrorPercent));

colormap(comparisonFigure, turbo);
comparisonColorbar = colorbar(improvedAxes);
comparisonColorbar.Label.String = '\epsilon_r';

function formatMapAxes(axesHandle, gridSize, colorLimits, titleText)
axis(axesHandle, 'equal');
xlim(axesHandle, [0.5 gridSize(1)+0.5]);
ylim(axesHandle, [0.5 gridSize(2)+0.5]);
set(axesHandle, 'YDir', 'normal');
clim(axesHandle, colorLimits);
xlabel(axesHandle, 'x index');
ylabel(axesHandle, 'y index');
title(axesHandle, titleText);
end

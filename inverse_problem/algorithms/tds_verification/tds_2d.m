% tds_2d Verify homogeneous target permittivity from transmission delays.
% Run build_cfg.m first to construct cfg in the current workspace.

assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'tds_2d:MissingConfig', 'Run build_cfg.m before tds_2d.m.');
assert(isfield(cfg, 'targets') && isscalar(cfg.targets), ...
    'tds_2d:ExpectedSingleTarget', ...
    'The verification requires exactly one configured target.');
assert(isfield(cfg.targets(1), 'mask') && ...
    isequal(size(cfg.targets(1).mask), [cfg.Nx cfg.Ny]), ...
    'tds_2d:InvalidTargetMask', ...
    'cfg.targets(1).mask must be a full-grid target mask.');

%% Verification controls
rayWidthPixels = 15;
maxAngle = deg2rad(40);
clampNegativeDelays = true;

targetMask = logical(cfg.targets(1).mask);
backgroundEpsr = cfg.grid.background.epsr;
knownTargetEpsr = cfg.targets(1).material.epsr;
numAntennas = cfg.antennas.numAntennas;

%% Pair geometry and the retained off-boresight filter
[distance, ~, betaTx, betaRx] = buildCircularPairGeometry( ...
    cfg.antennas.pos, cfg.dx, cfg.dy, cfg.antennas.center);
pairEligible = abs(betaTx) <= maxAngle & ...
    abs(betaRx) <= maxAngle & isfinite(distance);
pairEligible(eye(numAntennas) == 1) = false;

fprintf('Angle-eligible directed TX-RX pairs: %d out of %d.\n', ...
    nnz(pairEligible), numAntennas*(numAntennas - 1));

%% Object/reference simulations and excess-delay extraction
delaySteps = nan(numAntennas, numAntennas);
arrivalObjectSteps = nan(numAntennas, numAntennas);
arrivalReferenceSteps = nan(numAntennas, numAntennas);

simulationTimer = tic;
for tx = 1:numAntennas
    cfg.source.samples(:) = 0;
    cfg.source.samples(tx, :) = cfg.source.func((0:cfg.Nt-1) .* cfg.dt);

    rxList = setdiff(1:numAntennas, tx);
    receiverIndices = cfg.antennas.pos(rxList, :);

    objectResult = fdtd_mex(cfg);
    referenceCfg = fdtdmat.setBackgroundDefault(cfg, cfg.grid.background);
    referenceResult = fdtd_mex(referenceCfg);

    objectArrivals = getArrivalSteps( ...
        objectResult.Ez, receiverIndices, cfg);
    referenceArrivals = getArrivalSteps( ...
        referenceResult.Ez, receiverIndices, referenceCfg);

    arrivalObjectSteps(tx, rxList) = objectArrivals.';
    arrivalReferenceSteps(tx, rxList) = referenceArrivals.';
    delaySteps(tx, rxList) = objectArrivals.' - referenceArrivals.';

    clear objectResult referenceResult
    fprintf('Finished transmitter %d / %d.\n', tx, numAntennas);
end
simulationRunTime = toc(simulationTimer);

%% Delay-derived average permittivity for every directed pair
deltaTime = delaySteps .* cfg.dt;
if clampNegativeDelays
    deltaTime(deltaTime < 0) = 0;
end

averageEpsr = (sqrt(backgroundEpsr) + ...
    (cfg.c0 .* deltaTime) ./ distance).^2;
averageEpsr(eye(numAntennas) == 1) = NaN;

%% Proposed arithmetic pixel-mixture inversion
[rayEstimates, targetPixelCounts, backgroundPixelCounts, ...
    pairAccepted, rayCoverage] = evaluateTargetRays( ...
    averageEpsr, backgroundEpsr, cfg.antennas.pos, targetMask, ...
    rayWidthPixels, pairEligible);

validRayEstimates = rayEstimates(pairAccepted);
if isempty(validRayEstimates)
    error('tds_2d:NoValidTargetRays', ...
        ['No finite angle-eligible rays intersected the target. ' ...
        'Inspect the target position, arrivals, and ray width.']);
end

recoveredTargetEpsr = mean(validRayEstimates);
reconstruction = reconstructHomogeneousTargetMap( ...
    targetMask, backgroundEpsr, recoveredTargetEpsr);
groundTruth = reconstructHomogeneousTargetMap( ...
    targetMask, backgroundEpsr, knownTargetEpsr);

absoluteError = abs(recoveredTargetEpsr - knownTargetEpsr);
relativeErrorPercent = 100 * absoluteError / abs(knownTargetEpsr);

fprintf('\nTransmission-delay homogeneous-target verification\n');
fprintf('  Ray width:                    %d pixels\n', rayWidthPixels);
fprintf('  Valid directed target rays:   %d\n', nnz(pairAccepted));
fprintf('  Known target epsr:            %.6f\n', knownTargetEpsr);
fprintf('  Recovered target epsr:        %.6f\n', recoveredTargetEpsr);
fprintf('  Absolute error:               %.6f\n', absoluteError);
fprintf('  Relative error:               %.2f %%\n', relativeErrorPercent);
fprintf('  Per-ray epsr min/mean/std/max: %.6f / %.6f / %.6f / %.6f\n', ...
    min(validRayEstimates), mean(validRayEstimates), ...
    std(validRayEstimates), max(validRayEstimates));
fprintf('  FDTD simulation time:         %.2f s\n\n', simulationRunTime);

%% Accepted thick-ray coverage
coverageFigure = figure('Name', 'TD verification ray coverage', ...
    'NumberTitle', 'off');
coverageAxes = axes('Parent', coverageFigure);
imagesc(coverageAxes, 1:cfg.Nx, 1:cfg.Ny, rayCoverage.');
axis(coverageAxes, 'equal');
xlim(coverageAxes, [0.5 cfg.Nx + 0.5]);
ylim(coverageAxes, [0.5 cfg.Ny + 0.5]);
set(coverageAxes, 'YDir', 'normal');
xlabel(coverageAxes, 'x index');
ylabel(coverageAxes, 'y index');
title(coverageAxes, sprintf( ...
    'Accepted %d-pixel directed-ray coverage (%d rays)', ...
    rayWidthPixels, nnz(pairAccepted)));
colormap(coverageAxes, turbo);
coverageColorbar = colorbar(coverageAxes);
coverageColorbar.Label.String = 'Accepted directed-ray overlap count';
hold(coverageAxes, 'on');
contour(coverageAxes, 1:cfg.Nx, 1:cfg.Ny, double(targetMask.'), ...
    [0.5 0.5], 'w-', 'LineWidth', 2, 'DisplayName', 'Known target');
plot(coverageAxes, cfg.antennas.pos(:, 1), cfg.antennas.pos(:, 2), ...
    'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
    'DisplayName', 'Antennas');
legend(coverageAxes, 'Location', 'best');
hold(coverageAxes, 'off');

%% Ground-truth and reconstructed maps with a shared color scale
comparisonFigure = figure('Name', 'TD quantitative verification', ...
    'NumberTitle', 'off');
comparisonLayout = tiledlayout(comparisonFigure, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

colorLimits = [min([backgroundEpsr knownTargetEpsr recoveredTargetEpsr]), ...
    max([backgroundEpsr knownTargetEpsr recoveredTargetEpsr])];
if colorLimits(1) == colorLimits(2)
    colorLimits = colorLimits + [-0.5 0.5];
end

truthAxes = nexttile(comparisonLayout);
imagesc(truthAxes, 1:cfg.Nx, 1:cfg.Ny, groundTruth.');
formatMapAxes(truthAxes, cfg, colorLimits, ...
    sprintf('Ground truth, \\epsilon_r = %.4f', knownTargetEpsr));

reconstructionAxes = nexttile(comparisonLayout);
imagesc(reconstructionAxes, 1:cfg.Nx, 1:cfg.Ny, reconstruction.');
formatMapAxes(reconstructionAxes, cfg, colorLimits, ...
    sprintf('Recovered, \\epsilon_r = %.4f (error %.2f%%)', ...
    recoveredTargetEpsr, relativeErrorPercent));

colormap(comparisonFigure, turbo);
comparisonColorbar = colorbar(reconstructionAxes);
comparisonColorbar.Label.String = '\epsilon_r';

%% Local plotting helper
function formatMapAxes(axesHandle, cfg, colorLimits, titleText)
axis(axesHandle, 'equal');
xlim(axesHandle, [0.5 cfg.Nx + 0.5]);
ylim(axesHandle, [0.5 cfg.Ny + 0.5]);
set(axesHandle, 'YDir', 'normal');
clim(axesHandle, colorLimits);
xlabel(axesHandle, 'x index');
ylabel(axesHandle, 'y index');
title(axesHandle, titleText);
end

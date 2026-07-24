function measurement = acquireTransmissionMeasurements(cfg, options)
%acquireTransmissionMeasurements Run object/reference scans for one case.

if nargin < 2 || isempty(options)
    options = struct();
end

maximumAngleDegrees = getOption(options, 'MaximumAngleDegrees', 40);
maximumMatchedDelaySteps = getOption(options, ...
    'MaximumMatchedDelaySteps', 120);
windowHalfWidthSteps = getOption(options, 'WindowHalfWidthSteps', ...
    max(4, round(4*cfg.source.pulseWidth/cfg.dt)));

numAntennas = cfg.antennas.numAntennas;
[distance, ~, betaTx, betaRx] = buildCircularPairGeometry( ...
    cfg.antennas.pos, cfg.dx, cfg.dy, cfg.antennas.center);
maximumAngle = deg2rad(maximumAngleDegrees);
pairEligible = abs(betaTx) <= maximumAngle & ...
    abs(betaRx) <= maximumAngle & isfinite(distance);
pairEligible(eye(numAntennas) == 1) = false;

arrivalObjectSteps = nan(numAntennas);
arrivalReferenceSteps = nan(numAntennas);
prominentDelaySteps = nan(numAntennas);
matchedDelaySteps = nan(numAntennas);
matchedScores = nan(numAntennas);

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
    [matchedDelays, scores] = getMatchedEnvelopeDelaySteps( ...
        objectResult.Ez, referenceResult.Ez, receiverIndices, ...
        referenceCfg, windowHalfWidthSteps, maximumMatchedDelaySteps);

    arrivalObjectSteps(tx, rxList) = objectArrivals.';
    arrivalReferenceSteps(tx, rxList) = referenceArrivals.';
    prominentDelaySteps(tx, rxList) = ...
        objectArrivals.' - referenceArrivals.';
    matchedDelaySteps(tx, rxList) = matchedDelays.';
    matchedScores(tx, rxList) = scores.';

    clear objectResult referenceResult
    fprintf('Finished transmitter %d / %d.\n', tx, numAntennas);
end

measurement = struct();
measurement.version = 1;
measurement.simulationRunTime = toc(simulationTimer);
measurement.gridSize = [cfg.Nx cfg.Ny];
measurement.dx = cfg.dx;
measurement.dy = cfg.dy;
measurement.dt = cfg.dt;
measurement.c0 = cfg.c0;
measurement.backgroundEpsr = cfg.grid.background.epsr;
measurement.targetMask = logical(cfg.targets(1).mask);
measurement.targetEpsr = cfg.targets(1).material.epsr;
measurement.targetConductivity = cfg.targets(1).material.cond_e;
measurement.antennaPositions = cfg.antennas.pos;
measurement.antennaCenter = cfg.antennas.center;
measurement.distance = distance;
measurement.pairEligible = pairEligible;
measurement.maximumAngleDegrees = maximumAngleDegrees;
measurement.arrivalObjectSteps = arrivalObjectSteps;
measurement.arrivalReferenceSteps = arrivalReferenceSteps;
measurement.prominentDelaySteps = prominentDelaySteps;
measurement.matchedDelaySteps = matchedDelaySteps;
measurement.matchedScores = matchedScores;
measurement.windowHalfWidthSteps = windowHalfWidthSteps;
measurement.maximumMatchedDelaySteps = maximumMatchedDelaySteps;
end

function value = getOption(options, name, defaultValue)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

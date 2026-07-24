function scan = runFdtdDbimScan( ...
    cfg, frequencies, doiLinearIndex, includeDoiFields)
%runFdtdDbimScan Run one sequential circular-array FDTD scan.

if exist('fdtd_mex', 'file') ~= 3
    error('runFdtdDbimScan:MissingFdtdMex', ...
        'Build forward_solver/mex/fdtd_mex before running DBIM.');
end
if ~isstruct(cfg) || ~all(isfield(cfg, ...
        {'Nx', 'Ny', 'Nt', 'dt', 'snapshotStart', 'snapshotStride', ...
        'antennas', 'source'}))
    error('runFdtdDbimScan:InvalidConfig', ...
        'cfg is missing required FDTD scan fields.');
end
if ~isfield(cfg.antennas, 'pos') || ~isfield(cfg.antennas, 'numAntennas') || ...
        ~isfield(cfg.source, 'func') || ...
        ~isa(cfg.source.func, 'function_handle')
    error('runFdtdDbimScan:InvalidConfig', ...
        'cfg must contain antenna positions and source.func.');
end
if ~isscalar(includeDoiFields) || ...
        ~(islogical(includeDoiFields) || isnumeric(includeDoiFields))
    error('runFdtdDbimScan:InvalidIncludeDoiFields', ...
        'includeDoiFields must be a logical scalar.');
end
includeDoiFields = logical(includeDoiFields);
frequencies = double(frequencies(:).');
if isempty(frequencies) || any(~isfinite(frequencies)) || any(frequencies <= 0)
    error('runFdtdDbimScan:InvalidFrequencies', ...
        'frequencies must contain finite positive values.');
end

numAntennas = double(cfg.antennas.numAntennas);
if ~isequal(size(cfg.antennas.pos), [numAntennas 2])
    error('runFdtdDbimScan:InvalidAntennaPositions', ...
        'cfg.antennas.pos must be numAntennas-by-2.');
end
receiverTimes = (0:cfg.Nt-1) .* cfg.dt;
sourcePulse = double(cfg.source.func(receiverTimes));
sourcePulse = sourcePulse(:).';
if numel(sourcePulse) ~= cfg.Nt || ~isreal(sourcePulse) || ...
        any(~isfinite(sourcePulse))
    error('runFdtdDbimScan:InvalidSourceSamples', ...
        'cfg.source.func must return Nt finite real samples.');
end

doiLinearIndex = double(doiLinearIndex(:));
if includeDoiFields && (isempty(doiLinearIndex) || ...
        any(~isfinite(doiLinearIndex)) || ...
        any(doiLinearIndex ~= round(doiLinearIndex)) || ...
        any(doiLinearIndex < 1) || any(doiLinearIndex > cfg.Nx*cfg.Ny))
    error('runFdtdDbimScan:InvalidDoiIndex', ...
        'doiLinearIndex must contain valid grid indices.');
end

fourierRuntime = 0;
fourierTimer = tic;
sourceSpectrum = extractFrequencyComponents( ...
    sourcePulse, receiverTimes, frequencies);
fourierRuntime = fourierRuntime + toc(fourierTimer);
sourceMagnitude = abs(sourceSpectrum);
if max(sourceMagnitude) <= eps || ...
        any(sourceMagnitude < 1e-8 * max(sourceMagnitude))
    error('runFdtdDbimScan:InsufficientSourceSpectrum', ...
        'The source has insufficient energy at one or more requested frequencies.');
end

numFrequencies = numel(frequencies);
receiverFields = complex(zeros(numAntennas, numAntennas, numFrequencies));
if includeDoiFields
    doiFields = complex(zeros(numel(doiLinearIndex), ...
        numAntennas, numFrequencies));
else
    doiFields = complex(zeros(0, numAntennas, numFrequencies));
end

scanCfg = cfg;
scanCfg.returnRxSignals = true;
scanCfg.returnHx = false;
scanCfg.returnHy = false;
scanCfg.returnEz = includeDoiFields;
if includeDoiFields
    if scanCfg.snapshotStride <= 0
        error('runFdtdDbimScan:InvalidSnapshotStride', ...
            'A positive snapshotStride is required for DOI field output.');
    end
    firstSnapshot = max(scanCfg.snapshotStart, 0);
    snapshotSteps = firstSnapshot:scanCfg.snapshotStride:(scanCfg.Nt - 1);
    snapshotTimes = snapshotSteps .* scanCfg.dt;
else
    scanCfg.snapshotStride = 0;
    snapshotTimes = zeros(1, 0);
end
scanCfg.source.samples = zeros(numAntennas, cfg.Nt);

scanTimer = tic;
simulationRuntime = 0;
for transmitter = 1:numAntennas
    scanCfg.source.samples(:) = 0;
    transmitterPulse = double(scanCfg.source.func(receiverTimes));
    scanCfg.source.samples(transmitter, :) = transmitterPulse(:).';
    simulationTimer = tic;
    mexResult = fdtd_mex(scanCfg);
    simulationRuntime = simulationRuntime + toc(simulationTimer);

    if ~isfield(mexResult, 'rx_signals') || ...
            ~isequal(size(mexResult.rx_signals), [numAntennas cfg.Nt])
        error('runFdtdDbimScan:InvalidReceiverOutput', ...
            'fdtd_mex returned an unexpected receiver trace array.');
    end
    fourierTimer = tic;
    receiverSpectrum = extractFrequencyComponents( ...
        mexResult.rx_signals, receiverTimes, frequencies);
    receiverFields(:, transmitter, :) = ...
        reshape(receiverSpectrum, numAntennas, 1, numFrequencies);

    if includeDoiFields
        if ~isfield(mexResult, 'Ez') || isempty(mexResult.Ez)
            error('runFdtdDbimScan:MissingFieldOutput', ...
                'fdtd_mex did not return the requested Ez history.');
        end
        fieldHistory = reshape(mexResult.Ez, cfg.Nx*cfg.Ny, []);
        if size(fieldHistory, 2) ~= numel(snapshotTimes)
            error('runFdtdDbimScan:SnapshotCountMismatch', ...
                'The returned Ez history does not match the configured snapshots.');
        end
        doiHistory = fieldHistory(doiLinearIndex, :);
        doiSpectrum = extractFrequencyComponents( ...
            doiHistory, snapshotTimes, frequencies);
        doiFields(:, transmitter, :) = ...
            reshape(doiSpectrum, numel(doiLinearIndex), 1, numFrequencies);
    end
    fourierRuntime = fourierRuntime + toc(fourierTimer);

    clear mexResult fieldHistory doiHistory
    fprintf('DBIM FDTD scan: transmitter %d/%d complete.\n', ...
        transmitter, numAntennas);
end

scan = struct();
scan.receiverFields = receiverFields;
scan.doiFields = doiFields;
scan.sourceSpectrum = sourceSpectrum;
scan.frequencies = frequencies;
scan.snapshotTimes = snapshotTimes;
scan.runtime = toc(scanTimer);
scan.simulationRuntime = simulationRuntime;
scan.fourierRuntime = fourierRuntime;
end

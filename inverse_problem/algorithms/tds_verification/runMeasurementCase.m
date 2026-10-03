function measurement = runMeasurementCase(verificationSeed, forceRecompute)
%runMeasurementCase Acquire and cache one deterministic verification case.

if nargin < 1 || isempty(verificationSeed)
    verificationSeed = 1;
end
if nargin < 2 || isempty(forceRecompute)
    forceRecompute = false;
end
validateattributes(verificationSeed, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive', '<=', 2^32-1}, ...
    mfilename, 'verificationSeed');
validateattributes(forceRecompute, {'logical', 'numeric'}, ...
    {'scalar'}, mfilename, 'forceRecompute');

algorithmDirectory = fileparts(mfilename('fullpath'));
resultsDirectory = fullfile(algorithmDirectory, 'results');
if ~isfolder(resultsDirectory)
    mkdir(resultsDirectory);
end
cacheFile = fullfile(resultsDirectory, ...
    sprintf('measurement_seed_%04d.mat', verificationSeed));

if ~logical(forceRecompute) && isfile(cacheFile)
    cached = load(cacheFile, 'measurement');
    if isfield(cached, 'measurement') && ...
            isfield(cached.measurement, 'version') && ...
            cached.measurement.version == 1
        measurement = cached.measurement;
        fprintf('Loaded cached measurement for seed %d.\n', verificationSeed);
        return
    end
end

addpath(algorithmDirectory, '-begin');
setup = struct('targetOpts', struct(), 'randomSeed', verificationSeed);
cfg = build_cfg(setup);
measurement = acquireTransmissionMeasurements(cfg);
measurement.seed = verificationSeed;
save(cacheFile, 'measurement');
fprintf('Saved measurement for seed %d to %s.\n', ...
    verificationSeed, cacheFile);
end

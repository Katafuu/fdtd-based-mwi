function [results, cfg] = runFbts(cfg, options)
% Test double: one compute failure, otherwise deterministic solver timings.
assert(~isfile(fullfile(fileparts(options.outputDirectory), 'batch_timings.csv')), ...
    'fbtsTest:PrematureCsv', 'The CSV must not exist while cases are running.');
if isequal(cfg.antennas.originalIndices, 1)
    error('fbtsTest:ComputeFailure', 'Injected reconstruction failure.');
end
count = cfg.antennas.numAntennas;
perIteration = reshape(1:count*options.numIterations, count, options.numIterations);
fdtd = struct('measurement', (1:count)', 'forward', perIteration, ...
    'adjoint', perIteration + 1, 'sensitivity_baseline', perIteration + 2, ...
    'sensitivity_perturbation', perIteration + 3);
results = struct('total_runtime', 0, 'timing', struct('fdtd', fdtd));
end

% fbts_demo Reconstruct, plot, and save lossless relative permittivity by FBTS.
% Run build_cfg.m immediately before this script.
% Optional fbtsOptions: numIterations, sensitivityDownsampleFactor,
% outputDirectory, and runLabel. See resolveFbtsOptions for defaults.

assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'fbts:MissingConfig', 'Run build_cfg.m before fbts_demo.m.');
assert(exist('fdtd_mex', 'file') == 3, ...
    'fbts:MissingFdtdMex', ...
    'Build forward_solver/mex/fdtd_mex before running fbts_demo.m.');

if exist('fbtsOptions', 'var') == 1
    activeFbtsOptions = resolveFbtsOptions(fbtsOptions);
else
    activeFbtsOptions = resolveFbtsOptions();
end
[results, cfg] = runFbts(cfg, activeFbtsOptions);
[results, fbtsFigureHandles] = saveFbtsResults(cfg, results, activeFbtsOptions);

numIterations = results.num_iterations;
trueFigure = fbtsFigureHandles.trueFigure;
estimatedFigure = fbtsFigureHandles.estimatedFigure;
convergenceFigure = fbtsFigureHandles.convergenceFigure;
costFigure = fbtsFigureHandles.costFigure;

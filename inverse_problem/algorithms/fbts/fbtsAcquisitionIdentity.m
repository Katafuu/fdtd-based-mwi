function identity = fbtsAcquisitionIdentity(cfg)
%fbtsAcquisitionIdentity Inputs influencing a scan, independent of its subset.
[~, time, pulse, ~, source] = prepareFbtsSource(cfg);
physical = struct();
for name = {'Nx','Ny','sizeZ','dx','dy','dt','Nt','pml','init'}
    if isfield(cfg, name{1}), physical.(name{1}) = cfg.(name{1}); end
end
for name = {'epsr','cond_e','cond_m','murx','mury'}
    physical.grid.(name{1}) = cfg.grid.(name{1});
end
physical.source = source;
physical.time = time;
physical.pulse = pulse;
solverPath = which('fdtd_mex');
if isfield(cfg,'solverBackend') && strcmp(cfg.solverBackend,'cuda')
    solverPath = which('fdtd_cuda');
end
assert(~isempty(solverPath), 'fbts:MissingFdtdMex', 'Build fdtd_mex first.');
physical.solver_hash = fbtsHash(solverPath, 'file');
identity = struct('physics_hash', fbtsHash(physical), ...
    'solver_hash', physical.solver_hash, 'time', time, 'source', source);
end

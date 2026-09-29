function sample = timeFbtsCase(cfg,options,measurements)
%timeFbtsCase Worker-safe benchmark entry; compact result for screening.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(fileparts(folder)));
addpath(folder,fullfile(root,'buildLib'),fullfile(root,'forward_solver','mex'));
startedUtc=fbtsUtcNow(); timer=tic;
result=runFbts(cfg,options,measurements,struct('computeImageError',false));
sample=struct('started_utc',startedUtc,'finished_utc',fbtsUtcNow(), ...
    'wall_seconds',toc(timer),'timing',result.timing, ...
    'memory',fbtsMemorySnapshot(),'completed_iterations',result.num_iterations);
end

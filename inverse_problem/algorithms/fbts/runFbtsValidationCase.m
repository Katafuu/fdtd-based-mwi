function [result,states] = runFbtsValidationCase(cfg,options,measurements,directory)
%runFbtsValidationCase Preserve every optimizer state without changing solver math.
if ~isfolder(directory), mkdir(directory); end
checkpoint=fullfile(directory,'checkpoint.mat');
states=cell(1,options.numIterations);
execution=struct('computeImageError',false,'checkpointFile',checkpoint, ...
    'historyFile',fullfile(directory,'history_work.mat'),'onIteration',@capture);
result=runFbts(cfg,options,measurements,execution);
saveFbtsAtomic(fullfile(directory,'iteration_states.mat'),struct('states',{states}));
    function capture(iteration)
        saved=load(checkpoint,'checkpoint');
        states{iteration}=saved.checkpoint.state;
    end
end

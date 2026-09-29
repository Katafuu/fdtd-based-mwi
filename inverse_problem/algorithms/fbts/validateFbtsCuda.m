function report = validateFbtsCuda(prepared,directory)
%validateFbtsCuda Actual-device gate; never substitutes CPU for missing CUDA.
hardware = fbtsHardwareInfo();
assert(hardware.gpu_available && hardware.cuda_mex_available,'fbts:CudaUnavailable', ...
    'Run this gate on the NVIDIA PC after building fdtd_cuda.');
if ~isfolder(directory), mkdir(directory); end
report = struct('passed',false,'hardware',hardware,'cases',[]);
% Direct field/receiver comparison exercises both propagation signs before FBTS.
c=prepared.scenes{1}; c.Nt=512;
[c,~,pulse]=prepareFbtsSource(c);
c.source.samples=zeros(c.antennas.numAntennas,c.Nt); c.source.samples(1,:)=pulse;
c.returnEz=true; c.returnHx=false; c.returnHy=false; c.returnRxSignals=true;
c.snapshotStart=0; c.snapshotStride=1;
[x,y]=find(c.antennas.doiMask); c.ezRegion=[min(x),max(x),min(y),max(y)];
for reverse=[false true]
    c.solverBackend='cpu'; cpu=fbtsSolve(c,reverse);
    c.solverBackend='cuda'; gpu=fbtsSolve(c,reverse);
    for name={'Ez','rx_signals'}
        a=gpu.(name{1}); b=cpu.(name{1});
        assert(all(isfinite(a),'all') && norm(a(:)-b(:))<=1e-12+1e-8*norm(b(:)), ...
            'fbts:CudaMismatch','CUDA direct field/receiver comparison failed.');
    end
    saveFbtsAtomic(fullfile(directory,sprintf('fields_reverse_%d.mat',reverse)),struct('cpu',cpu,'cuda',gpu));
end
clear cpu gpu
ids = [1 1 1 4 7 10]; sets = {1:8,[2 4 6 8],6,6,6,6};
for k = 1:numel(ids)
    cfg = subsetFbtsConfig(prepared.scenes{ids(k)},sets{k});
    options = struct('numIterations',3,'sensitivityDownsampleFactor',prepared.sensitivity.factor);
    [reference,referenceStates] = runFbtsValidationCase(cfg,options,[], ...
        fullfile(directory,sprintf('control_%02d_cpu',k)));
    options.solverBackend = 'cuda';
    [actual,actualStates] = runFbtsValidationCase(cfg,options,[], ...
        fullfile(directory,sprintf('control_%02d_cuda',k)));
    entry = struct('scene_index',ids(k),'active_indices',sets{k},'errors',struct());
    groups = {{'Ez_measured','rx_model_history'}, ...
        {'raw_gradient_epsr','gradient_epsr','rx_sensitivity_coarse_history','Ez_sensitivity', ...
        'alpha','step_a','step_q','pr_beta'}, ...
        {'epsr_history_doi','cost'}};
    tolerance = [1e-8 1e-7 1e-6];
    for group = 1:3
        for name = groups{group}
            a=actual.(name{1}); b=reference.(name{1});
            difference=norm(a(:)-b(:)); scale=norm(b(:));
            entry.errors.(name{1})=difference/max(scale,realmin);
            assert(all(isfinite(a),'all') && difference <= 1e-12+tolerance(group)*scale, ...
                'fbts:CudaMismatch','CUDA comparison failed for %s in control %d.',name{1},k);
        end
    end
    assert(isequal(actual.pr_restart,reference.pr_restart) && ...
        isequal(actual.projection_applied,reference.projection_applied), ...
        'fbts:CudaMismatch','CUDA changed restart/projection decisions.');
    for iteration=1:options.numIterations
        for name={'gradEpsr','projectedGradEpsr','directionEpsr'}
            a=actualStates{iteration}.(name{1}); b=referenceStates{iteration}.(name{1});
            assert(all(isfinite(a),'all') && norm(a(:)-b(:))<=1e-12+1e-7*norm(b(:)), ...
                'fbts:CudaMismatch','CUDA gradient/direction trajectory differs.');
        end
    end
    report.cases = [report.cases;entry]; %#ok<AGROW>
    saveFbtsAtomic(fullfile(directory,sprintf('control_%02d.mat',k)), ...
        struct('reference',reference,'actual',actual,'comparison',entry));
end
report.passed=true;
saveFbtsAtomic(fullfile(directory,'cuda_validation.mat'),struct('report',report));
end

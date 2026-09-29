function result = fbtsSolve(cfg,tr)
%fbtsSolve Explicit backend dispatch; never silently changes a running study.
if nargin < 2, tr = false; end
backend = 'cpu';
if isfield(cfg,'solverBackend'), backend = cfg.solverBackend; end
switch backend
    case 'cpu'
        result = fdtd_mex(cfg,tr);
    case 'cuda'
        assert(exist('fdtd_cuda','file')==3 && exist('gpuDevice','file')==2, ...
            'fbts:CudaUnavailable','Build fdtd_cuda and install Parallel Computing Toolbox on the NVIDIA PC.');
        index = 1; if isfield(cfg,'gpuDeviceIndex'), index = cfg.gpuDeviceIndex; end
        device = gpuDevice;
        if device.Index ~= index, device = gpuDevice(index); end
        result = fdtd_cuda(cfg,tr);
        wait(device); % Returned solve time includes device completion.
    otherwise
        error('fbts:InvalidBackend','Unknown solver backend.');
end
end

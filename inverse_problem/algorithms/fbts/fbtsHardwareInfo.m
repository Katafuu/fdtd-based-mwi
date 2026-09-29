function info = fbtsHardwareInfo()
%fbtsHardwareInfo Read-only preflight; suitable for sharing with benchmark data.
info = struct('matlab',version,'platform',computer,'parallel_available',fbtsParallelAvailable(), ...
    'logical_processors',java.lang.Runtime.getRuntime().availableProcessors(), ...
    'cpu_mex_available',exist('fdtd_mex','file')==3,'cuda_mex_available',exist('fdtd_cuda','file')==3, ...
    'gpu_available',false,'gpu_message','Parallel Computing Toolbox GPU functions are not installed.');
compiler = mex.getCompilerConfigurations('C','Selected');
if isempty(compiler), info.c_compiler = 'Not configured'; else, info.c_compiler = compiler.Name; end
compiler = mex.getCompilerConfigurations('C++','Selected');
if isempty(compiler), info.cpp_compiler = 'Not configured'; else, info.cpp_compiler = compiler.Name; end
info.mexcuda_available=exist('mexcuda','file')==2;
if isfile('/proc/meminfo')
    content = fileread('/proc/meminfo');
    for key = {'MemTotal','MemAvailable'}
        value = regexp(content,[key{1} ':\s+(\d+) kB'],'tokens','once');
        if ~isempty(value), info.([key{1} '_bytes']) = 1024*str2double(value{1}); end
    end
elseif ispc
    [~,system] = memory; info.MemTotal_bytes = system.PhysicalMemory.Total;
    info.MemAvailable_bytes = system.PhysicalMemory.Available;
end
info.cpu_model = 'unknown';
if isfile('/proc/cpuinfo')
    cpu = fileread('/proc/cpuinfo');
    model = regexp(cpu,'model name\s*:\s*([^\n]+)','tokens','once');
    if ~isempty(model), info.cpu_model = strtrim(model{1}); end
    packages = regexp(cpu,'physical id\s*:\s*(\d+)','tokens');
    cores = regexp(cpu,'core id\s*:\s*(\d+)','tokens');
    if ~isempty(cores) && numel(packages)==numel(cores)
        ids = cellfun(@(p,c) [p{1} ':' c{1}],packages,cores,'UniformOutput',false);
        info.physical_cores = min(numel(unique(ids)),info.logical_processors);
    end
elseif ispc
    info.cpu_model = getenv('PROCESSOR_IDENTIFIER');
end
% Conservative fallback when the physical topology is unavailable.
if ~isfield(info,'physical_cores'), info.physical_cores = max(1,floor(info.logical_processors/2)); end
if exist('gpuDevice','file')==2
    try
        d = gpuDevice;
        info.gpu_available = true; info.gpu_message = '';
        info.gpu = struct('name',d.Name,'index',d.Index,'compute_capability',d.ComputeCapability, ...
            'total_memory_bytes',d.TotalMemory,'available_memory_bytes',d.AvailableMemory);
        for name={'DriverVersion','ToolkitVersion'}
            if isprop(d,name{1}), info.gpu.(name{1})=d.(name{1}); end
        end
    catch exception
        info.gpu_message = exception.message;
    end
end
identity={computer,getenv('HOSTNAME'),getenv('COMPUTERNAME'),info.cpu_model, ...
    info.logical_processors,info.physical_cores};
if isfield(info,'gpu'), identity{end+1}=rmfield(info.gpu,'available_memory_bytes'); end
info.machine_id=fbtsHash(identity);
if nargout == 0, disp(info); end
end

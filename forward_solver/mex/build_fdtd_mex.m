% build_fdtd_mex Build the TMz FDTD MATLAB MEX interface.
%
% The MEX gateway accepts canonical MATLAB [x,y] material, PML, antenna,
% source, and initial-state arrays. Returned fields use the same canonical
% orientation; internal C storage conversion is handled by the gateway.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);

sources = {
    fullfile(repoRoot, 'mex', 'fdtd_mex.c')
    fullfile(repoRoot, 'fdtd', 'solver.c')
    fullfile(repoRoot, 'fdtd', 'solver_tr.c')
    fullfile(repoRoot, 'fdtd', 'gridtmz.c')
    fullfile(repoRoot, 'fdtd', 'updatetmz.c')
    fullfile(repoRoot, 'fdtd', 'updatetmz_tr.c')
    fullfile(repoRoot, 'fdtd', 'pml.c')
    fullfile(repoRoot, 'fdtd', 'pml_tr.c')
    fullfile(repoRoot, 'fdtd', 'ricker.c')
    fullfile(repoRoot, 'fdtd', 'abctmz.c')
    fullfile(repoRoot, 'fdtd', 'grid1dez.c')
    fullfile(repoRoot, 'fdtd', 'tfsftmz.c')
    fullfile(repoRoot, 'utility', 'fdtd-config.c')
    fullfile(repoRoot, 'utility', 'fdtd-io.c')
    fullfile(repoRoot, 'utility', 'snapshot2d.c')
};

includeArgs = {
    ['-I' fullfile(repoRoot, 'utility')]
    ['-I' fullfile(repoRoot, 'fdtd')]
};

outDir = scriptDir;
% Optional workspace value: struct('mode','openmp','outputDirectory',...).
buildMode = 'serial';
if exist('fdtdBuildOptions','var')
    if isfield(fdtdBuildOptions,'mode'), buildMode = char(fdtdBuildOptions.mode); end
    if isfield(fdtdBuildOptions,'outputDirectory'), outDir = char(fdtdBuildOptions.outputDirectory); end
end
buildMode = validatestring(buildMode,{'serial','openmp','cuda'});
flags = {};
if strcmp(buildMode,'openmp')
    if ispc
        selected = mex.getCompilerConfigurations('C','Selected');
        if contains(selected.Name,'Microsoft')
            flags = {'COMPFLAGS=$COMPFLAGS /openmp'};
        else
            flags = {'CFLAGS=$CFLAGS -fopenmp','LDFLAGS=$LDFLAGS -fopenmp'};
        end
    elseif isunix && ~ismac
        flags = {'CFLAGS=$CFLAGS -fopenmp','LDFLAGS=$LDFLAGS -fopenmp'};
    else
        error('fdtd:UnsupportedOpenMPBuild','OpenMP build is configured for Linux GCC and Windows compilers.');
    end
end
if ~isfolder(outDir), mkdir(outDir); end
if strcmp(buildMode,'cuda')
    assert(exist('mexcuda','file')==2,'fdtd:CudaToolboxMissing','Install Parallel Computing Toolbox on the NVIDIA PC.');
    clear fdtd_cuda
    flags = {'-DFDTD_CUDA','NVCCFLAGS=$NVCCFLAGS --fmad=false -arch=sm_75'};
    mexcuda('-R2018a',flags{:},includeArgs{:},sources{:}, ...
        fullfile(repoRoot,'fdtd','cuda_solver.cu'),'-outdir',outDir,'-output','fdtd_cuda');
    compiler = mex.getCompilerConfigurations('C++','Selected');
    infoFile = 'fdtd_cuda_build_info.mat';
else
    clear fdtd_mex
    mex('-R2018a', flags{:}, includeArgs{:}, sources{:}, '-outdir', outDir, '-output', 'fdtd_mex');
    compiler = mex.getCompilerConfigurations('C','Selected');
    infoFile = 'fdtd_build_info.mat';
end
buildInfo = struct('mode',buildMode,'compiler',compiler.Name,'compiler_version',compiler.Version, ...
    'flags',{flags},'matlab_version',version,'platform',computer);
save(fullfile(outDir,infoFile),'buildInfo');

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
mex('-R2018a', includeArgs{:}, sources{:}, '-outdir', outDir, '-output', 'fdtd_mex');

%[text] # Run TMz FDTD Through Files
%[text] This Live Script runs the same shared case through the standalone executable. It writes the input CSV/config files, launches |tmzdemo2.exe|, and reads the output CSV back into MATLAB.
%[text] Expected workspace outputs are |outputPath|, |snapshots|, |lastEzFromFile|, and |standaloneRuntime|.
exampleDir = fileparts(mfilename('fullpath'));
if isempty(exampleDir)
    exampleDir = pwd;
end
solverRoot = fileparts(exampleDir);
workspaceRoot = fileparts(solverRoot);
addpath(fullfile(workspaceRoot, 'buildLib'));
%%
%[text] ## Load or Build the Case
%[text] The file-backed runner uses the normal MATLAB-orientation |grid| and |pml| arrays created by |build_fdtd_case|. If those variables are missing, this section builds the default case first.
if ~exist('cfg', 'var') || ~isstruct(cfg) || ...
        ~exist('grid', 'var') || ~exist('pml', 'var') || ...
        ~isfield(cfg, 'Nx') || ~isfield(cfg, 'Ny') || ~isfield(cfg, 'Nt') || ...
        ~isfield(cfg, 'antennas') || ~isfield(cfg, 'source')
    run(fullfile(exampleDir, 'build_fdtd_case.m'));
end
%%
%[text] ## Prepare Output Folders
%[text] The standalone executable reads scalar and matrix input files from |results/fdtd_input| and writes the Ez snapshot CSV into |results|.
resultsDir = fullfile(solverRoot, 'results');
inputDir = fullfile(resultsDir, 'fdtd_input');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
if ~exist(inputDir, 'dir')
    mkdir(inputDir);
end
%%
%[text] ## Write Scalar Configuration
%[text] The config file stores scalar settings and the antenna count. Antenna coordinates remain one-based, matching MATLAB.
fileSnapshotStart = 0;
fileSnapshotStride = 1;
configPath = fullfile(inputDir, 'grid_config.txt');
configLines = [
    "Nx " + cfg.Nx
    "Ny " + cfg.Ny
    "Nt " + cfg.Nt
    "dx " + compose('%.17g', cfg.dx)
    "dy " + compose('%.17g', cfg.dy)
    "dt " + compose('%.17g', cfg.dt)
    "num_antennas " + cfg.antennas.numAntennas
    "ax " + compose('%.17g', cfg.pml.ax)
    "ay " + compose('%.17g', cfg.pml.ay)
    "az " + compose('%.17g', cfg.pml.az)
    "pml_type " + string(cfg.pml.type)
    "pml_enabled " + double(cfg.pml.enabled)
    "snapshot_start " + fileSnapshotStart
    "snapshot_stride " + fileSnapshotStride
];
writelines(configLines, configPath);
%%
%[text] ## Write Antennas and Sources
%[text] The standalone files use the same public matrix shapes as cfg.
assert(isequal(size(cfg.antennas.pos), [cfg.antennas.numAntennas 2]), ...
    'cfg.antennas.pos must be numAntennas-by-2.');
assert(isequal(size(cfg.source.samples), [cfg.antennas.numAntennas cfg.Nt]), ...
    'cfg.source.samples must be numAntennas-by-Nt.');
writematrix(cfg.antennas.pos, fullfile(inputDir, 'antennas.csv'));
writematrix(cfg.source.samples, fullfile(inputDir, 'source.csv'));
%%
%[text] ## Write Material and PML Maps
%[text] These CSV files preserve the standalone row-wise file contract while reading canonical maps directly from |cfg|.
writematrix(cfg.grid.epsr, fullfile(inputDir, 'epsr.csv'));
writematrix(cfg.grid.murx, fullfile(inputDir, 'murx.csv'));
writematrix(cfg.grid.mury, fullfile(inputDir, 'mury.csv'));
writematrix(cfg.grid.cond_e, fullfile(inputDir, 'cond_e.csv'));
writematrix(cfg.grid.cond_m, fullfile(inputDir, 'cond_m.csv'));
writematrix(cfg.pml.condx, fullfile(inputDir, 'condx.csv'));
writematrix(cfg.pml.condy, fullfile(inputDir, 'condy.csv'));
writematrix(cfg.pml.kx, fullfile(inputDir, 'kx.csv'));
writematrix(cfg.pml.ky, fullfile(inputDir, 'ky.csv'));
%%
%[text] ## Locate or Build the Standalone Executable
%[text] If |tmzdemo2.exe| is missing, this section invokes the existing Makefile before running the file-backed simulation.
exePath = fullfile(solverRoot, 'tmzdemo2.exe');
if ~isfile(exePath)
    oldFolder = pwd;
    cd(solverRoot);
    [buildStatus, buildLog] = system('mingw32-make');
    cd(oldFolder);
    assert(buildStatus == 0, buildLog);
end
%%
%[text] ## Run the Standalone Solver
%[text] The executable receives the input directory and output CSV path explicitly. A zero process status means the output file is ready to read.
outputPath = fullfile(resultsDir, 'ez_field.csv');
runCommand = sprintf('"%s" "%s" "%s"', exePath, inputDir, outputPath);
tic;
[runStatus, runLog] = system(runCommand);
standaloneRuntime = toc;
assert(runStatus == 0, runLog);
%%
%[text] ## Read Output Snapshots
%[text] The first output row stores grid dimensions. Remaining rows are flattened Ez snapshots, which are reshaped into |lastEzFromFile| for inspection by a separate visualization script.
outputLines = readlines(outputPath);
gridSize = sscanf(outputLines(1), '%d,%d').';
snapshots = readmatrix(outputPath, 'NumHeaderLines', 1);
expectedSnapshotColumns = prod(gridSize);
assert(size(snapshots, 2) == expectedSnapshotColumns, ...
    'Snapshot row width does not match Nx*Ny from the output header.');
lastEzFromFile = reshape(snapshots(end, :), [gridSize(2), gridSize(1)]).';
%%
%[text] ## Workspace Result
%[text] The file-backed run is complete. The raw snapshot matrix and final Ez field are available for separate visualization or comparison against the MEX output.
disp("Standalone solver finished. snapshots and lastEzFromFile are available in the workspace.");
fprintf("Standalone runtime: %.6f seconds\n", standaloneRuntime);
%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
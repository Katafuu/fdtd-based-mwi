%[text] # Run TMz FDTD Through MEX
%[text] This Live Script runs the shared case through the MATLAB MEX interface. It keeps the simulation data in MATLAB memory and leaves |Ez| in the workspace for separate visualization.
%[text] Expected workspace outputs are |fdtdMexResult|, |Ez|, |rxSignals|, and |mexRuntime|.
exampleDir = fileparts(mfilename('fullpath'));
if isempty(exampleDir)
    exampleDir = pwd;
end
solverRoot = fileparts(exampleDir);
workspaceRoot = fileparts(solverRoot);
addpath(fullfile(workspaceRoot, 'buildLib')); %[output:1d5063b3]
addpath(fullfile(solverRoot, 'mex')); %[output:43faf56e]
%%
%[text] ## Load or Build the Case
%[text] The MEX runner depends on the shared |cfg| struct created by |build\_fdtd\_case|. If |cfg| is not already in the workspace, this section builds the default case first.
if ~exist('cfg', 'var') || ~isstruct(cfg) || ...
        ~isfield(cfg, 'grid') || ~isfield(cfg, 'pml') || ...
        ~isfield(cfg, 'antennas') || ~isfield(cfg, 'source') || ...
        ~isfield(cfg.source, 'samples')
    run(fullfile(exampleDir, 'build_fdtd_case.m'));
end
%%
%[text] ## Locate or Build the MEX Function
%[text] The compiled MEX file lives in the |mex| folder. If MATLAB cannot find |fdtd\_mex|, this section runs the MEX build script before launching the simulation.
if exist('fdtd_mex', 'file') ~= 3
    run(fullfile(solverRoot, 'mex', 'build_fdtd_mex.m'));
end
%%
%[text] ## Run the Solver
%[text] The MEX call passes |cfg| directly into C. The C solver uses the same update loop as the standalone executable, but output is captured in MATLAB arrays instead of written to CSV.
tic;
fdtdMexResult = fdtd_mex(cfg); %[output:4a356fe9]
mexRuntime = toc;
%%
%[text] ## Expose Returned Ez
%[text] The MEX gateway returns |Ez| directly in canonical |Nx|-by-|Ny| MATLAB orientation.
if isfield(fdtdMexResult, 'Ez') && ~isempty(fdtdMexResult.Ez)
    Ez = fdtdMexResult.Ez;
else
    Ez = [];
end
%%
%[text] ## Expose Receiver Signals
%[text] Receiver traces already use the public |numAntennas|-by-|Nt| layout and therefore require no transposition.
if isfield(fdtdMexResult, 'rx_signals') && ...
        ~isempty(fdtdMexResult.rx_signals)
    rxSignals = fdtdMexResult.rx_signals;
else
    rxSignals = [];
end
%%
%[text] ## Workspace Result
%[text] The final field and optional receiver traces are now available as |Ez| and |rxSignals|. Any plotting, comparison, or post-processing should happen in a separate visualization script.
disp("MEX solver finished. Ez and rxSignals are available in the workspace.");
fprintf("MEX runtime: %.6f seconds\n", mexRuntime);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:1d5063b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Name is nonexistent or not a directory: \/buildLib"}}
%---
%[output:43faf56e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Name is nonexistent or not a directory: \/tmp\/mex"}}
%---
%[output:4a356fe9]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'fdtd_mex'."}}
%---

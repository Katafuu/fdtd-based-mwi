function [tr_result, focusPoint] = tr_mwi(cfg)
% tr_mwi Run time reversal and return its focus location in grid indices.

%% Load relevant libraries

algorithmDir = fileparts(mfilename('fullpath'));
algorithmsDir = fileparts(algorithmDir);
inverseDir = fileparts(algorithmsDir);
workspaceRoot = fileparts(inverseDir);

buildLibDir = fullfile(workspaceRoot, 'buildLib');
trLibDir = fullfile(algorithmsDir, 'time_reversal', 'lib');
mwiLibDir = fullfile(algorithmsDir, 'transmission_mwi', 'lib');
hybridLibDir = fullfile(algorithmDir, 'lib');
mexDir = fullfile(workspaceRoot, 'forward_solver', 'mex');

addpath(buildLibDir, '-end');
addpath(trLibDir, '-end');
addpath(mwiLibDir, '-end');
addpath(hybridLibDir, '-begin');
addpath(mexDir, '-end');

tr_result = itr_run(cfg);
[~, focusLinearIndex] = max(abs(tr_result.focusMagFrame(:))); [focusX, focusY] = ind2sub(size(tr_result.focusMagFrame), focusLinearIndex);
focusPoint = [focusX, focusY];

end

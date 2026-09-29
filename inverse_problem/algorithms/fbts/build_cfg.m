% build_cfg Build one FBTS scene, optionally using fbtsBuildOptions.
% Example: fbtsBuildOptions = struct('seed',42,'numAntennas',12);
scriptDirectory = fileparts(mfilename('fullpath'));
addpath(scriptDirectory);
if exist('fbtsBuildOptions','var')
    cfg = buildFbtsConfig(fbtsBuildOptions);
else
    cfg = buildFbtsConfig();
end

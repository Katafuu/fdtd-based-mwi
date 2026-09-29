function prepared = prepareFbtsStudy(directory,solverThreads)
%prepareFbtsStudy Validate physics and persist the 12-scene study before running.
if nargin < 2, solverThreads = 1; end
if ~isfolder(directory), mkdir(directory); end
preparedFile = fullfile(directory,'prepared.mat'); provenance = fbtsCodeProvenance();
if isfile(preparedFile)
    loaded = load(preparedFile,'prepared'); prepared = loaded.prepared;
    assert(strcmp(prepared.code_hash,provenance.code_hash),'fbts:StudyMismatch', ...
        'Code/build changed after study preparation. Use a new study directory.');
    return
end
timer = tic;
window = calibrateFbtsStudyWindow(fullfile(directory,'calibration','window'),solverThreads);
sceneFile = fullfile(directory,'scenes.mat');
if isfile(sceneFile)
    saved = load(sceneFile);
    assert(strcmp(saved.code_hash,provenance.code_hash),'fbts:StudyMismatch','Saved scene code/build differs.');
    scenes = saved.scenes; index = saved.index; generationTiming = saved.timing;
else
    [scenes,index,generationTiming] = buildFbtsStudyScenes(window.recording_seconds,42);
    % Accepted coordinates and RNG states are authoritative on resume.
    saveFbtsAtomic(sceneFile,struct('scenes',{scenes},'index',index,'timing',generationTiming, ...
        'code_hash',provenance.code_hash));
end
sensitivity = calibrateFbtsSensitivity(scenes,fullfile(directory,'calibration','sensitivity'),solverThreads);
prepared = struct('schema_version',1,'scenes',{scenes},'scene_index',index, ...
    'window',window,'sensitivity',sensitivity,'generation_timing',generationTiming, ...
    'code_hash',provenance.code_hash,'hardware',fbtsHardwareInfo(), ...
    'num_iterations',15,'selection_seed',123,'configurations_per_scene',29, ...
    'total_configurations',348,'normal_fdtd_calls',86496,'wall_seconds',toc(timer));
saveFbtsAtomic(preparedFile,struct('prepared',prepared));
writetable(index,fullfile(directory,'scenes.csv'));
end

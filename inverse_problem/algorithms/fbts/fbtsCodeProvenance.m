function provenance = fbtsCodeProvenance(destination)
%fbtsCodeProvenance Content identity and optional source snapshot, including untracked code.
fbts = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(fbts)));
folders = {fbts, fullfile(root,'buildLib'), fullfile(root,'forward_solver')};
files = {};
for k = 1:numel(folders)
    listing = dir(fullfile(folders{k},'**','*'));
    for j = 1:numel(listing)
        if listing(j).isdir, continue; end
        [~,~,ext] = fileparts(listing(j).name);
        if ~ismember(ext,{'.m','.c','.h','.cu','.cuh'}) && ~strcmp(listing(j).name,'Makefile'), continue; end
        path = fullfile(listing(j).folder,listing(j).name);
        if any(contains(path, { ...
                [filesep 'tests' filesep],[filesep 'figs' filesep], ...
                [filesep 'results' filesep],[filesep 'provenance' filesep]})), continue; end
        files{end+1} = path; %#ok<AGROW>
    end
end
files = sort(files);
entries = repmat(struct('path','','sha256',''),numel(files),1);
for k = 1:numel(files)
    relative = files{k}(numel(root)+2:end);
    entries(k) = struct('path',relative,'sha256',fbtsHash(files{k},'file'));
    if nargin > 0 && ~isempty(destination)
        target = fullfile(destination,relative);
        if ~isfolder(fileparts(target)), mkdir(fileparts(target)); end
        copyfile(files{k},target);
    end
end
provenance = struct('files',entries,'source_hash',fbtsHash(entries), ...
    'solver_hash',fbtsHash(which('fdtd_mex'),'file'), ...
    'matlab_version',version,'platform',computer, ...
    'compiler','unknown (existing MEX)','build_flags','unknown (existing MEX)');
provenance.cuda_solver_hash = '';
if exist('fdtd_cuda','file')==3, provenance.cuda_solver_hash = fbtsHash(which('fdtd_cuda'),'file'); end
buildFile = fullfile(fileparts(which('fdtd_mex')),'fdtd_build_info.mat');
if isfile(buildFile)
    build = load(buildFile,'buildInfo'); provenance.cpu_build = build.buildInfo;
end
if exist('fdtd_cuda','file')==3
    buildFile = fullfile(fileparts(which('fdtd_cuda')),'fdtd_cuda_build_info.mat');
    if isfile(buildFile)
        build = load(buildFile,'buildInfo'); provenance.cuda_build = build.buildInfo;
    end
end
provenance.code_hash = fbtsHash({provenance.source_hash,provenance.solver_hash,provenance.cuda_solver_hash});
% Query only repository facts, never the environment or credentials.
old = pwd; restore = onCleanup(@() cd(old)); %#ok<NASGU>
cd(root);
[ok,revision] = system('git rev-parse HEAD');
if ok == 0, provenance.git_revision = strtrim(revision); else, provenance.git_revision = 'unknown'; end
[ok,status] = system('git status --porcelain');
if ok == 0, provenance.git_status = status; else, provenance.git_status = 'unknown'; end
end

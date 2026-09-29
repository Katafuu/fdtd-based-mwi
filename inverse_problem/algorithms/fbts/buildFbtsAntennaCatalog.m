function catalog = buildFbtsAntennaCatalog(cfg, symmetry, disabledCounts, seed, sceneId, maxCandidates)
%buildFbtsAntennaCatalog One reproducible uniform representative per orbit.
if nargin < 2, symmetry = 'dihedral'; end
n = cfg.antennas.numAntennas;
if nargin < 3 || isempty(disabledCounts), disabledCounts = 0:n-1; end
if nargin < 4, seed = 0; end
if nargin < 5, sceneId = 'scene_000001'; end
if nargin < 6, maxCandidates = 2e6; end
symmetry = validatestring(symmetry, {'dihedral','rotation','none'});
validateattributes(seed, {'numeric'}, {'scalar','integer','nonnegative','<=',2^32-1});
validateattributes(disabledCounts, {'numeric'}, {'vector','integer','>=',0,'<',n});
if ~strcmp(symmetry, 'none'), validateRing(cfg); end
counts = unique(disabledCounts(:).');
candidateCount = sum(arrayfun(@(k) nchoosek(n,k), counts));
assert(candidateCount <= maxCandidates, 'fbts:CatalogTooLarge', ...
    'Enumerating %.0f candidates exceeds maxCandidates=%.0f. Narrow disabledCounts or raise the limit explicitly.', ...
    candidateCount, maxCandidates);
template = struct('class_id','','num_antennas',n,'disabled_count',0, ...
    'canonical_gaps',[],'canonical_disabled',[],'orbit_size',0, ...
    'active_indices',[],'disabled_indices',[],'rotation',0,'reflected',false, ...
    'selection_seed',0,'selection_rng_state',[],'symmetry',symmetry);
catalog = repmat(template, 0, 1);
seen = containers.Map('KeyType','char','ValueType','logical');
for k = counts
    disabled = 1:k;
    while true
        [key,gaps] = fbtsAntennaClassKey(n, disabled, symmetry);
        if ~isKey(seen, key)
            seen(key) = true;
            [orbit, transforms] = orbitOf(n,disabled,symmetry);
            h = fbtsHash({double(seed),char(sceneId),key});
            localSeed = hex2dec(h(1:8));
            stream = RandStream('mt19937ar','Seed',localSeed);
            before = stream.State;
            choice = randi(stream,size(orbit,1));
            entry = template;
            entry.class_id = key; entry.disabled_count = k;
            entry.canonical_gaps = gaps; entry.canonical_disabled = disabled;
            entry.orbit_size = size(orbit,1);
            entry.disabled_indices = orbit(choice,:);
            entry.active_indices = setdiff(1:n,entry.disabled_indices);
            entry.rotation = transforms(choice,1); entry.reflected = logical(transforms(choice,2));
            entry.selection_seed = localSeed; entry.selection_rng_state = before;
            catalog(end+1,1) = entry; %#ok<AGROW>
        end
        % Stream lexicographic combinations instead of constructing nchoosek arrays.
        j = find(disabled < n-k+(1:k),1,'last');
        if isempty(j), break; end
        disabled(j) = disabled(j)+1;
        disabled(j+1:k) = disabled(j)+(1:k-j);
    end
end
end

function [orbit, transforms] = orbitOf(n,d,symmetry)
if isempty(d) || strcmp(symmetry,'none')
    orbit = reshape(d,1,[]); transforms = [0 0]; return
end
signs = 1;
if strcmp(symmetry,'dihedral'), signs = [1 -1]; end
orbit = zeros(n*numel(signs),numel(d));
transforms = zeros(size(orbit,1),2); row = 0;
for s = signs
    for r = 0:n-1
        row = row+1;
        orbit(row,:) = sort(mod(s*(d-1)+r,n)+1);
        transforms(row,:) = [r,s==-1];
    end
end
[orbit,first] = unique(orbit,'rows','sorted');
transforms = transforms(first,:);
end

function validateRing(cfg)
a = cfg.antennas; n = a.numAntennas;
assert(isfield(a,'arrayType') && strcmp(a.arrayType,'circular') && ...
    isfield(a,'center') && isfield(a,'radius') && cfg.dx == cfg.dy, ...
    'fbts:NonuniformArray', 'Symmetry pruning requires a uniformly spaced circular array with dx=dy. Use symmetry="none" otherwise.');
theta = (0:n-1)'*2*pi/n;
expected = round(a.center + a.radius*[cos(theta),sin(theta)]);
assert(isequal(a.pos,expected) && size(unique(a.pos,'rows'),1)==n, ...
    'fbts:NonuniformArray', 'Antenna positions/order do not match the nominal circular array.');
end

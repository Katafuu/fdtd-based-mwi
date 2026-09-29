function [key, gaps] = fbtsAntennaClassKey(n, disabled, symmetry)
%fbtsAntennaClassKey Ordered circular gaps modulo rotations and reflections.
if nargin < 3, symmetry = 'dihedral'; end
symmetry = validatestring(symmetry, {'dihedral','rotation','none'});
validateattributes(n, {'numeric'}, {'scalar','integer','positive','finite'});
d = sort(disabled(:).');
assert(all(d >= 1 & d <= n & d == fix(d)) && numel(unique(d)) == numel(d), ...
    'fbts:InvalidSubset', 'Disabled indices must be distinct original antenna IDs.');
k = numel(d);
if k == 0
    gaps = [];
else
    gaps = diff([d, n+d(1)]);
    candidates = zeros(k, k);
    for j = 1:k, candidates(j,:) = circshift(gaps, [0,j-1]); end
    if strcmp(symmetry, 'dihedral')
        candidates = [candidates; fliplr(candidates)];
    end
    candidates = sortrows(candidates);
    gaps = candidates(1,:);
end
values = gaps;
if strcmp(symmetry, 'none'), values = d; end
width = numel(num2str(n));
key = sprintf('n%d_d%0*d_%s', n, width, k, symmetry);
for j = 1:numel(values), key = [key sprintf('_%0*d',width,values(j))]; end %#ok<AGROW>
end

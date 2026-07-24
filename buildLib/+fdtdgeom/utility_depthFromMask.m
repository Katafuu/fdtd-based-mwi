function depth = utility_depthFromMask(mask, maxDepth, metric)
%utility_depthFromMask Compute distance inside a logical mask.

if nargin < 2 || isempty(maxDepth)
    maxDepth = inf;
end
if nargin < 3 || isempty(metric)
    metric = "euclidean";
end
if ~islogical(mask)
    mask = mask ~= 0;
end
if ~isscalar(maxDepth) || ~isnumeric(maxDepth) || maxDepth <= 0
    error('fdtdgeom:utility_depthFromMask:InvalidMaxDepth', 'maxDepth must be a positive scalar.');
end

metric = char(string(metric));
if exist('bwdist', 'file') == 2
    depth = bwdist(~mask, metric);
else
    depth = bruteForceDistance(mask, metric);
end

depth(~mask) = 0;
depth = min(depth, maxDepth);
end

function depth = bruteForceDistance(mask, metric)
depth = zeros(size(mask));
[outsideRows, outsideCols] = find(~mask);
[insideRows, insideCols] = find(mask);

if isempty(insideRows)
    return;
end
if isempty(outsideRows)
    depth(mask) = inf;
    return;
end

for pointIndex = 1:numel(insideRows)
    dr = insideRows(pointIndex) - outsideRows;
    dc = insideCols(pointIndex) - outsideCols;
    switch lower(metric)
        case 'euclidean'
            distance = hypot(dr, dc);
        case 'cityblock'
            distance = abs(dr) + abs(dc);
        case 'chessboard'
            distance = max(abs(dr), abs(dc));
        case 'quasi-euclidean'
            distance = hypot(dr, dc);
        otherwise
            error('fdtdgeom:utility_depthFromMask:InvalidMetric', ...
                'metric must be euclidean, cityblock, chessboard, or quasi-euclidean.');
    end
    depth(insideRows(pointIndex), insideCols(pointIndex)) = min(distance);
end
end

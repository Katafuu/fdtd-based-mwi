function h = plotXY(A, xRange, yRange)
% plotXY Plot an internally stored A(x,y) matrix with x horizontal.
%
% MATLAB imagesc displays matrix rows vertically and columns horizontally.
% If A is stored as A(x,y), transpose for display only.

    if nargin < 2 || isempty(xRange)
        xRange = 1:size(A, 1);
    end
    if nargin < 3 || isempty(yRange)
        yRange = 1:size(A, 2);
    end

    h = imagesc(xRange, yRange, A.');
    axis equal tight;
    set(gca, 'YDir', 'normal');
    xlabel('x index');
    ylabel('y index');
    colorbar;
end

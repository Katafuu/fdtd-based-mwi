function newgrid = BlurAndDownsample(highGrid, refine)
%BlurAndDownsample FFT-blur high-resolution material maps and sample to a coarse grid.

if nargin ~= 2
    error('fdtdmat:BlurAndDownsample:InvalidInputCount', ...
        'Expected highGrid and refine inputs.');
end
if ~isstruct(highGrid)
    error('fdtdmat:BlurAndDownsample:InvalidGrid', 'highGrid must be a grid struct.');
end
if ~isnumeric(refine) || ~isscalar(refine) || ~isfinite(refine) || refine < 1 || fix(refine) ~= refine
    error('fdtdmat:BlurAndDownsample:InvalidRefine', ...
        'refine must be a positive integer scalar.');
end
if ~isfield(highGrid, 'sizeXY') || ~isfield(highGrid, 'spacingXY') || ...
        ~isfield(highGrid, 'originPhysical') || ~isfield(highGrid, 'background')
    error('fdtdmat:BlurAndDownsample:InvalidGrid', ...
        'highGrid must contain sizeXY, spacingXY, originPhysical, and background fields.');
end

coarseSize = highGrid.sizeXY ./ refine;
if any(mod(highGrid.sizeXY, refine) ~= 0)
    error('fdtdmat:BlurAndDownsample:InvalidRefine', ...
        'highGrid.sizeXY must be divisible by refine.');
end

newgrid = fdtdmat.createGrid( ...
    coarseSize, ...
    highGrid.spacingXY .* refine, ...
    highGrid.background, ...
    highGrid.originPhysical);

fieldsToSmooth = {'epsr', 'cond_e', 'cond_m'};
for fieldIndex = 1:numel(fieldsToSmooth)
    fieldName = fieldsToSmooth{fieldIndex};
    if isfield(highGrid, fieldName)
        newgrid.(fieldName) = blurAndSample(highGrid.(fieldName), refine, coarseSize, fieldName);
    end
end
end

function coarseMap = blurAndSample(highMap, refine, coarseSize, fieldName)
if ~isnumeric(highMap) || ~ismatrix(highMap)
    error('fdtdmat:BlurAndDownsample:InvalidMap', ...
        'highGrid.%s must be a numeric 2-D matrix.', fieldName);
end

kernel = zeros(size(highMap));
kernel(1:refine, 1:refine) = 1 / refine^2;
kernel = circshift(kernel, -floor([refine refine] / 2));

convolved = real(ifft2(fft2(highMap) .* fft2(kernel)));

sampleX = floor(refine / 2) + 1 : refine : size(highMap, 1);
sampleY = floor(refine / 2) + 1 : refine : size(highMap, 2);
coarseMap = convolved(sampleX, sampleY);

if ~isequal(size(coarseMap), coarseSize)
    error('fdtdmat:BlurAndDownsample:DownsampleSizeMismatch', ...
        'Downsampled %s map does not match the expected coarse size.', fieldName);
end
end

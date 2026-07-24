function components = extractFrequencyComponents(samples, sampleTimes, frequencies)
%extractFrequencyComponents Evaluate selected Fourier components.
%   The final dimension of samples is interpreted as time. The transform
%   convention is sum(E(t).*exp(-1i*omega*t))*dt.

if ~isnumeric(samples) || isempty(samples) || any(~isfinite(samples(:)))
    error('extractFrequencyComponents:InvalidSamples', ...
        'samples must be a nonempty finite numeric array.');
end

sampleTimes = double(sampleTimes(:));
frequencies = double(frequencies(:).');
if any(~isfinite(sampleTimes)) || numel(sampleTimes) < 2
    error('extractFrequencyComponents:InvalidSampleTimes', ...
        'sampleTimes must contain at least two finite values.');
end
if isempty(frequencies) || any(~isfinite(frequencies)) || any(frequencies <= 0)
    error('extractFrequencyComponents:InvalidFrequencies', ...
        'frequencies must contain finite positive values.');
end

sampleSize = size(samples);
numTimeSamples = sampleSize(end);
if numel(sampleTimes) ~= numTimeSamples
    error('extractFrequencyComponents:TimeSizeMismatch', ...
        'The final samples dimension must match sampleTimes.');
end

timeSteps = diff(sampleTimes);
sampleInterval = mean(timeSteps);
uniformTolerance = 1e-10 * max(1, abs(sampleInterval));
if sampleInterval <= 0 || any(abs(timeSteps - sampleInterval) > uniformTolerance)
    error('extractFrequencyComponents:NonuniformSampleTimes', ...
        'sampleTimes must be strictly increasing and uniformly spaced.');
end

fourierKernel = exp(-1i * 2*pi * sampleTimes * frequencies);
sampleMatrix = reshape(samples, [], numTimeSamples);
componentMatrix = (sampleMatrix * fourierKernel) * sampleInterval;

outputSize = [sampleSize(1:end-1), numel(frequencies)];
components = reshape(componentMatrix, outputSize);
end

function info = fbtsMemorySnapshot()
%fbtsMemorySnapshot Process memory, not the sum of overlapping worker durations.
info = struct('resident_bytes',NaN,'peak_resident_bytes',NaN);
if isfile('/proc/self/status')
    value = fileread('/proc/self/status');
    current = regexp(value,'VmRSS:\s+(\d+) kB','tokens','once');
    peak = regexp(value,'VmHWM:\s+(\d+) kB','tokens','once');
    if ~isempty(current), info.resident_bytes=1024*str2double(current{1}); end
    if ~isempty(peak), info.peak_resident_bytes=1024*str2double(peak{1}); end
elseif ispc
    try
        process=System.Diagnostics.Process.GetCurrentProcess();
        info.resident_bytes=double(process.WorkingSet64);
        info.peak_resident_bytes=double(process.PeakWorkingSet64);
    catch
        % Unknown memory is explicitly NaN, never zero.
    end
end
end

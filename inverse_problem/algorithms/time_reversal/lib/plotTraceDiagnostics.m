function fig = plotTraceDiagnostics(traceInfo, dt)
%plotTraceDiagnostics Plot raw, windowed, and reversed traces separately.

if nargin < 2 || isempty(dt)
    dt = traceInfo.time(2) - traceInfo.time(1);
end
timeNs = (0:size(traceInfo.raw, 2)-1) .* dt .* 1e9;
fig = gobjects(1, 3);

fig(1) = figure('Name', 'Raw scattered receiver traces', 'Position', [100 100 1000 420]);
plot(timeNs, traceInfo.raw.')
grid on
xlabel('Time [ns]')
ylabel('E_z scattered')
title('Raw scattered receiver traces')

fig(2) = figure('Name', 'Windowed TR input traces', 'Position', [140 140 1000 420]);
plot(timeNs, traceInfo.windowed.')
grid on
xlabel('Time [ns]')
ylabel('E_z windowed')
title('Windowed TR input traces')

fig(3) = figure('Name', 'Time-reversed injection traces', 'Position', [180 180 1000 420]);
plot(timeNs, traceInfo.reversed.')
grid on
xlabel('Time [ns]')
ylabel('E_z^{TR}')
title('Time-reversed injection traces')
end
function fig = plotTrComparison(cfg, tr_result, plotOpts)
%plotTrComparison Plot target map and TR focus images in separate figures.

if nargin < 3 || isempty(plotOpts)
    plotOpts = struct();
end

txAntenna = getTxAntenna(cfg, plotOpts);
focusMax = max([tr_result.focus_mag_image(:); tr_result.focus_entropy_image(:)]);
fig = gobjects(1, 3);

fig(1) = figure('Name', 'Actual target map', 'Position', [100 100 760 650]);
imagesc(cfg.grid.epsr')
axis equal tight
set(gca, 'YDir', 'normal')
colorbar
title('Actual target map: \epsilon_r')
xlabel('x grid index')
ylabel('y grid index')
hold on
plotDomainOverlays(gca, cfg, txAntenna, 'k');
hold off

fig(2) = figure('Name', 'Max magnitude TR focus', 'Position', [140 120 760 650]);
imagesc(tr_result.focus_mag_image')
axis equal tight
set(gca, 'YDir', 'normal')
colorbar
if focusMax > 0 && isfinite(focusMax)
    clim([0 focusMax])
end
title(sprintf('Max magnitude focus, step %d', tr_result.focus_mag_step))
xlabel('x grid index')
ylabel('y grid index')
hold on
plotDomainOverlays(gca, cfg, txAntenna, 'w');
hold off

fig(3) = figure('Name', 'Minimum R TR focus', 'Position', [180 140 760 650]);
imagesc(tr_result.focus_entropy_image')
axis equal tight
set(gca, 'YDir', 'normal')
colorbar
if focusMax > 0 && isfinite(focusMax)
    clim([0 focusMax])
end
title(sprintf('Minimum R focus, step %d', tr_result.focus_entropy_step))
xlabel('x grid index')
ylabel('y grid index')
hold on
plotDomainOverlays(gca, cfg, txAntenna, 'w');
hold off
end

function txAntenna = getTxAntenna(cfg, plotOpts)
txAntenna = getOption(plotOpts, 'txAntenna', []);
if isempty(txAntenna) && isfield(cfg, 'antennas') && ...
        isfield(cfg.antennas, 'txAntennas') && ~isempty(cfg.antennas.txAntennas)
    txAntenna = cfg.antennas.txAntennas(1);
end
if isempty(txAntenna)
    txAntenna = 1;
end
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end
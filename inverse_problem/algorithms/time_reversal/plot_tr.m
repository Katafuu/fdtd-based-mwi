% plot_tr Plot time-reversal results already present in the workspace.

assert(exist('cfg', 'var') == 1 && isstruct(cfg), ...
    'plot_tr:MissingConfig', 'Run itr_simple.m before plot_tr.m.');
assert(exist('tr_result', 'var') == 1 && isstruct(tr_result), ...
    'plot_tr:MissingResult', 'Run itr_simple.m before plot_tr.m.');

trFigures = plotTrComparison(cfg, tr_result);

% run_experiment_suite Reproduce the eight-case estimator comparison.

verificationSeeds = 1:8;
forceRecompute = false;
measurements = cell(size(verificationSeeds));
for caseIndex = 1:numel(verificationSeeds)
    measurements{caseIndex} = runMeasurementCase( ...
        verificationSeeds(caseIndex), forceRecompute);
end

experimentResults = analyzeVerificationCases(measurements);
disp(experimentResults);

developmentRows = experimentResults.seed <= 5;
holdoutRows = experimentResults.seed >= 6;
fprintf('\nDevelopment seeds 1-5\n');
fprintf('  Original MAPE:       %.3f%%\n', ...
    mean(experimentResults.baselineErrorPercent(developmentRows)));
fprintf('  Maximum-chord MAPE:  %.3f%%\n', ...
    mean(experimentResults.maximumChordErrorPercent(developmentRows)));
fprintf('  Maximum-chord worst: %.3f%%\n', ...
    max(experimentResults.maximumChordErrorPercent(developmentRows)));
fprintf('\nUnseen holdout seeds 6-8\n');
fprintf('  Original MAPE:       %.3f%%\n', ...
    mean(experimentResults.baselineErrorPercent(holdoutRows)));
fprintf('  Maximum-chord MAPE:  %.3f%%\n', ...
    mean(experimentResults.maximumChordErrorPercent(holdoutRows)));
fprintf('  Maximum-chord worst: %.3f%%\n', ...
    max(experimentResults.maximumChordErrorPercent(holdoutRows)));
fprintf('\nAll eight cases\n');
fprintf('  Original MAPE:       %.3f%%\n', ...
    mean(experimentResults.baselineErrorPercent));
fprintf('  Maximum-chord MAPE:  %.3f%%\n', ...
    mean(experimentResults.maximumChordErrorPercent));
fprintf('  Maximum-chord worst: %.3f%%\n\n', ...
    max(experimentResults.maximumChordErrorPercent));

errorFigure = figure('Name', 'Estimator error across seeds', ...
    'NumberTitle', 'off');
errorAxes = axes('Parent', errorFigure);
bar(errorAxes, experimentResults.seed, ...
    [experimentResults.baselineErrorPercent ...
    experimentResults.maximumChordErrorPercent], 'grouped');
xlabel(errorAxes, 'Verification seed');
ylabel(errorAxes, 'Absolute relative error (%)');
title(errorAxes, 'Original versus maximum-chord reconstruction error');
legend(errorAxes, {'Original 15-pixel arithmetic', ...
    'Matched-delay maximum chord'}, 'Location', 'northwest');
grid(errorAxes, 'on');

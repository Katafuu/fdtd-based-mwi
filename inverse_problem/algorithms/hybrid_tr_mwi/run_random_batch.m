function summary = run_random_batch(numCases, firstSeed, targetOpts)
%run_random_batch Run seeded random hybrid reconstructions and summarize DOI error.
% Example: summary = run_random_batch(10, 1);
% A new results/batch_XXXX folder holds the LaTeX report and batch MAT files.
% Each tds_2d case retains its own figs/tds_2d/run_XXXX folder.

if nargin < 1 || isempty(numCases), numCases = 10; end
if nargin < 2 || isempty(firstSeed), firstSeed = 1; end
if nargin < 3 || isempty(targetOpts)
    targetOpts = struct( ...
        'numTargetsRange', [1 1], ...
        'allowedShapes', ["circle" "rectangle" "triangle"], ...
        'radiusRange', [10 28], ...
        'sideRange', [18 50]);
end
validateattributes(numCases, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'});
validateattributes(firstSeed, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'});
assert(firstSeed + numCases - 1 <= 2^32 - 1, ...
    'run_random_batch:SeedRange', 'Case seeds exceed the supported range.');
assert(isstruct(targetOpts) && isscalar(targetOpts), ...
    'run_random_batch:TargetOpts', 'targetOpts must be a scalar struct.');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir, '-begin');
resultsRoot = fullfile(scriptDir, 'results');
if ~isfolder(resultsRoot), mkdir(resultsRoot); end
existingBatches = dir(fullfile(resultsRoot, 'batch_*'));
batchNumbers = zeros(0, 1);
for idx = 1:numel(existingBatches)
    if ~existingBatches(idx).isdir, continue; end
    token = regexp(existingBatches(idx).name, ...
        '^batch_(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        batchNumbers(end + 1, 1) = str2double(token{1}); %#ok<AGROW>
    end
end
batchDir = fullfile(resultsRoot, ...
    sprintf('batch_%04d', max([0; batchNumbers]) + 1));
mkdir(batchDir);

cases = struct('caseNumber', {}, 'seed', {}, 'runDir', {}, ...
    'reconstructionFile', {}, 'parametersFile', {});
failures = struct('caseNumber', {}, 'seed', {}, ...
    'identifier', {}, 'message', {});
for caseNumber = 1:numCases
    seed = firstSeed + caseNumber - 1;
    fprintf('\nBatch case %d/%d, random seed %d.\n', ...
        caseNumber, numCases, seed);
    try
        cfg = build_cfg(struct( ...
            'targetOpts', targetOpts, 'randomSeed', seed));
        caseParameters = collectParameters(cfg, seed);
        result = tds_2d(cfg);
        runDir = result.outputDir;
        parametersFile = fullfile(runDir, 'case_parameters.mat');
        save(parametersFile, 'caseParameters');
        cases(end + 1) = struct( ... %#ok<AGROW>
            'caseNumber', caseNumber, 'seed', seed, ...
            'runDir', runDir, ...
            'reconstructionFile', fullfile(runDir, 'reconstruction.mat'), ...
            'parametersFile', parametersFile);
    catch caseError
        failures(end + 1) = struct( ... %#ok<AGROW>
            'caseNumber', caseNumber, 'seed', seed, ...
            'identifier', caseError.identifier, ...
            'message', caseError.message);
        warning('run_random_batch:CaseFailed', ...
            'Case %d (seed %d) failed: %s', ...
            caseNumber, seed, caseError.message);
    end
end

manifestFile = fullfile(batchDir, 'batch_manifest.mat');
save(manifestFile, 'cases', 'failures', 'numCases', ...
    'firstSeed', 'targetOpts');

% Read saved MAT files after all simulations; no in-memory result is used
% to calculate the batch scores.
scores = struct('caseNumber', {}, 'seed', {}, 'runName', {}, ...
    'runDir', {}, 'parameters', {}, 'estimatedEpsr', {}, ...
    'selectedChords', {}, 'mae', {}, 'rmse', {}, 'mape', {});
for idx = 1:numel(cases)
    savedReconstruction = load(cases(idx).reconstructionFile, 'result');
    savedParameters = load(cases(idx).parametersFile, 'caseParameters');
    savedResult = savedReconstruction.result;
    parameters = savedParameters.caseParameters;
    doiMask = logical(parameters.doiMask);
    truth = savedResult.trueEpsr;
    reconstruction = savedResult.reconstructedEpsr;
    assert(isequal(size(doiMask), size(truth), size(reconstruction)) ...
        && any(doiMask(:)), 'run_random_batch:InvalidSavedDoi', ...
        'Saved DOI, truth, and reconstruction must share a nonempty grid.');
    assert(all(truth(doiMask) > 0), ...
        'run_random_batch:InvalidTruth', ...
        'True epsr must be positive throughout the DOI.');
    difference = reconstruction(doiMask) - truth(doiMask);
    scores(end + 1) = struct( ... %#ok<AGROW>
        'caseNumber', cases(idx).caseNumber, ...
        'seed', cases(idx).seed, ...
        'runName', string(localRunName(cases(idx).runDir)), ...
        'runDir', cases(idx).runDir, ...
        'parameters', parameters, ...
        'estimatedEpsr', savedResult.recoveredTargetEpsr, ...
        'selectedChords', ...
            savedResult.maximumChordDetails.numSelectedDirectedRays, ...
        'mae', mean(abs(difference)), ...
        'rmse', sqrt(mean(difference.^2)), ...
        'mape', 100 * mean(abs(difference) ./ truth(doiMask)));
end

summary = struct('batchDir', batchDir, 'requestedCases', numCases, ...
    'successfulCases', numel(scores), 'failedCases', numel(failures), ...
    'scores', scores, 'failures', failures, ...
    'mae', metricSummary([scores.mae], scores), ...
    'rmse', metricSummary([scores.rmse], scores), ...
    'mape', metricSummary([scores.mape], scores));
save(fullfile(batchDir, 'batch_summary.mat'), 'summary');
writeLatexReport(fullfile(batchDir, 'batch_report.tex'), summary, targetOpts);

fprintf('\nSaved %d/%d successful cases and batch report in %s.\n', ...
    numel(scores), numCases, batchDir);
if ~isempty(scores)
    fprintf('DOI MAE: mean %.6f, min %.6f, max %.6f.\n', ...
        summary.mae.mean, summary.mae.minimum, summary.mae.maximum);
end
end

function parameters = collectParameters(cfg, seed)
targets = repmat(struct('name', "", 'properties', struct(), ...
    'epsr', NaN, 'condE', NaN), 1, numel(cfg.targets));
for idx = 1:numel(cfg.targets)
    targets(idx).name = string(cfg.targets(idx).name);
    targets(idx).properties = cfg.targets(idx).properties;
    targets(idx).epsr = cfg.targets(idx).material.epsr;
    targets(idx).condE = cfg.targets(idx).material.cond_e;
end
parameters = struct( ...
    'seed', seed, 'gridSize', [cfg.Nx cfg.Ny], ...
    'dx', cfg.dx, 'dy', cfg.dy, 'dt', cfg.dt, 'Nt', cfg.Nt, ...
    'numAntennas', cfg.antennas.numAntennas, ...
    'antennaRadiusCells', cfg.antennas.radius, ...
    'backgroundEpsr', cfg.grid.background.epsr, ...
    'backgroundCondE', cfg.grid.background.cond_e, ...
    'pulseWidth', cfg.source.pulseWidth, ...
    'doiMask', logical(cfg.antennas.doiMask), ...
    'targets', targets);
end

function name = localRunName(runDir)
[~, name] = fileparts(runDir);
end

function stats = metricSummary(values, scores)
if isempty(values)
    stats = struct('mean', NaN, 'minimum', NaN, ...
        'minimumRun', "", 'maximum', NaN, 'maximumRun', "");
    return;
end
[minimum, minIdx] = min(values);
[maximum, maxIdx] = max(values);
stats = struct('mean', mean(values), ...
    'minimum', minimum, 'minimumRun', scores(minIdx).runName, ...
    'maximum', maximum, 'maximumRun', scores(maxIdx).runName);
end

function writeLatexReport(reportFile, summary, targetOpts)
fid = fopen(reportFile, 'w');
assert(fid ~= -1, 'run_random_batch:ReportOpen', ...
    'Could not create %s.', reportFile);
fileCleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '\\documentclass[11pt]{article}\n');
fprintf(fid, '\\usepackage[margin=0.8in]{geometry}\n');
fprintf(fid, '\\usepackage{amsmath,booktabs,longtable}\n');
fprintf(fid, '\\begin{document}\n');
fprintf(fid, '\\section*{Hybrid TR--transmission-delay random-case results}\n');
fprintf(fid, 'Requested cases: %d. Successful: %d. Failed: %d.\\par\n', ...
    summary.requestedCases, summary.successfulCases, summary.failedCases);
fprintf(fid, 'Each case uses seeded random target geometry and estimated shape support. ');
fprintf(fid, 'DOI errors compare the saved reconstructed and true relative-permittivity ');
fprintf(fid, 'matrices only inside the saved DOI mask. Failed cases are excluded from ');
fprintf(fid, 'the statistics below; each successful case has equal weight.\\par\n');
fprintf(fid, '\\[\\mathrm{MAE}=\\operatorname{mean}_{p\\in\\mathrm{DOI}}');
fprintf(fid, '|\\hat\\epsilon_r(p)-\\epsilon_r(p)|,\\quad ');
fprintf(fid, '\\mathrm{MAPE}=100\\operatorname{mean}_{p\\in\\mathrm{DOI}}');
fprintf(fid, '\\frac{|\\hat\\epsilon_r(p)-\\epsilon_r(p)|}{\\epsilon_r(p)}.\\]\n');

if ~isempty(summary.scores)
    p = summary.scores(1).parameters;
    fprintf(fid, '\\subsection*{Shared system model}\n');
    fprintf(fid, '\\begin{tabular}{ll}\\toprule\n');
    fprintf(fid, 'Grid & $%d\\times%d$ cells, $\\Delta x=%.3g$ m, ', ...
        p.gridSize(1), p.gridSize(2), p.dx);
    fprintf(fid, '$\\Delta y=%.3g$ m\\\\\n', p.dy);
    fprintf(fid, 'Time & $\\Delta t=%.4f$ ps, $N_t=%d$\\\\\n', ...
        1e12 * p.dt, p.Nt);
    fprintf(fid, 'Array & %d antennas, radius %.0f cells\\\\\n', ...
        p.numAntennas, p.antennaRadiusCells);
    fprintf(fid, 'Background & $\\epsilon_r=%.4g$, ', p.backgroundEpsr);
    fprintf(fid, '$\\sigma=%.4g$ S/m\\\\\n', p.backgroundCondE);
    fprintf(fid, 'DOI & %d cells\\\\\n', nnz(p.doiMask));
    fprintf(fid, 'Source pulse width & %.4f ps\\\\\n', ...
        1e12 * p.pulseWidth);
    fprintf(fid, '\\bottomrule\\end{tabular}\\par\n');
end
fprintf(fid, 'Random-generation options: %s target(s); shapes %s; ', ...
    optionRange(targetOpts, 'numTargetsRange', [1 1]), ...
    strjoin(string(optionValue(targetOpts, 'allowedShapes', ...
        ["rectangle"])), ', '));
fprintf(fid, 'radius %s cells; side %s cells; target ', ...
    optionRange(targetOpts, 'radiusRange', [10 28]), ...
    optionRange(targetOpts, 'sideRange', [31 31]));
fprintf(fid, '$\\epsilon_r$ %s; $\\sigma$ %s S/m.\\par\n', ...
    optionRange(targetOpts, 'epsrRange', [2.5 6]), ...
    optionRange(targetOpts, 'condRange', [0.02 0.15]));

if ~isempty(summary.scores)
    fprintf(fid, '\\subsection*{DOI error summary}\n');
    fprintf(fid, '\\begin{tabular}{lrrrrr}\\toprule\n');
    fprintf(fid, 'Metric & Mean & Min & Min run & Max & Max run\\\\\\midrule\n');
    writeMetricRow(fid, 'MAE ($\epsilon_r$)', summary.mae);
    writeMetricRow(fid, 'RMSE ($\epsilon_r$)', summary.rmse);
    writeMetricRow(fid, 'MAPE (\%)', summary.mape);
    fprintf(fid, '\\bottomrule\\end{tabular}\n');

    fprintf(fid, '\\subsection*{Per-case system model}\n');
    fprintf(fid, '\\small\\begin{longtable}{rrllp{5.0cm}rr}\n');
    fprintf(fid, '\\toprule Case & Seed & Run & Shape & Geometry (grid cells) & ');
    fprintf(fid, 'True $\\epsilon_r$ & $\\sigma$ (S/m)\\\\\\midrule\\endhead\n');
    for idx = 1:numel(summary.scores)
        score = summary.scores(idx);
        for targetIdx = 1:numel(score.parameters.targets)
            target = score.parameters.targets(targetIdx);
            fprintf(fid, '%d & %d & %s & %s & %s & %.4f & %.4f\\\\\n', ...
                score.caseNumber, score.seed, ...
                latexUnderscores(score.runName), target.name, ...
                targetGeometry(target), target.epsr, target.condE);
        end
    end
    fprintf(fid, '\\bottomrule\\end{longtable}\\normalsize\n');

    fprintf(fid, '\\subsection*{Per-case reconstruction}\n');
    fprintf(fid, '\\begin{longtable}{lrrrrr}\n');
    fprintf(fid, '\\toprule Run & Estimated $\\epsilon_r$ & Chords & ');
    fprintf(fid, 'DOI MAE & DOI RMSE & DOI MAPE (\\%%)\\\\\\midrule\\endhead\n');
    for idx = 1:numel(summary.scores)
        score = summary.scores(idx);
        fprintf(fid, '%s & %.4f & %d & %.6f & %.6f & %.4f\\\\\n', ...
            latexUnderscores(score.runName), score.estimatedEpsr, ...
            score.selectedChords, score.mae, score.rmse, score.mape);
    end
    fprintf(fid, '\\bottomrule\\end{longtable}\n');
    fprintf(fid, 'Each run has its own figures, DOI text file, and MAT data under ');
    fprintf(fid, '\\texttt{figs/tds\\_2d/run\\_XXXX/} ');
    fprintf(fid, 'relative to the hybrid\\_tr\\_mwi folder.\\par\n');
end
if ~isempty(summary.failures)
    fprintf(fid, '\\subsection*{Failed cases}\n');
    fprintf(fid, '\\begin{tabular}{rrl}\\toprule Case & Seed & Error ID\\\\\\midrule\n');
    for idx = 1:numel(summary.failures)
        failure = summary.failures(idx);
        fprintf(fid, '%d & %d & %s\\\\\n', failure.caseNumber, ...
            failure.seed, latexUnderscores(failure.identifier));
    end
    fprintf(fid, '\\bottomrule\\end{tabular}\n');
end
fprintf(fid, '\\end{document}\n');
end

function writeMetricRow(fid, label, stats)
fprintf(fid, '%s & %.6f & %.6f & %s & %.6f & %s\\\\\n', ...
    label, stats.mean, stats.minimum, latexUnderscores(stats.minimumRun), ...
    stats.maximum, latexUnderscores(stats.maximumRun));
end

function value = optionValue(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end

function value = optionRange(opts, name, defaultValue)
range = optionValue(opts, name, defaultValue);
value = sprintf('[%.4g, %.4g]', range(1), range(2));
end

function value = latexUnderscores(value)
value = strrep(char(string(value)), '_', '\_');
end

function value = targetGeometry(target)
properties = target.properties;
switch string(target.name)
    case "circle"
        value = sprintf('center (%.0f, %.0f), radius %.0f', ...
            properties.center(1), properties.center(2), ...
            properties.radius);
    case "rectangle"
        bounds = properties.bounds;
        value = sprintf('x [%g, %g], y [%g, %g]', bounds);
    case "triangle"
        vertices = properties.vertices;
        value = sprintf('(%.1f, %.1f), (%.1f, %.1f), (%.1f, %.1f)', ...
            vertices(1, 1), vertices(1, 2), ...
            vertices(2, 1), vertices(2, 2), ...
            vertices(3, 1), vertices(3, 2));
    otherwise
        value = 'unknown';
end
end

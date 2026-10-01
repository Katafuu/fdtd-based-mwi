classdef testHybridTargetConfiguration < matlab.unittest.TestCase
    %testHybridTargetConfiguration Cover target setup and material variation.

    methods (Test)
        function fixedDefaultAndSeededRandomTargets(testCase)
            restoreState = configureHybridBuilder(); %#ok<NASGU>
            rng(1);
            fixedA = build_cfg();
            rng(999);
            fixedB = build_cfg(struct());
            testCase.verifyEqual(fixedA.targets, fixedB.targets);
            testCase.verifyEqual(fixedA.grid.epsr, fixedB.grid.epsr);
            testCase.verifyEqual(fixedA.pml, fixedB.pml);
            testCase.verifyEqual(fixedA.antennas, fixedB.antennas);
            testCase.verifyEqual(fixedA.targets(1).properties.center, [215 205]);
            testCase.verifyEqual(fixedA.targets(1).properties.radius, 20);
            testCase.verifyEqual(fixedA.targets(1).material.epsr, 5);
            testCase.verifyEqual(fixedA.targets(1).material.cond_e, 0.1);

            setup = struct('targetOpts', struct( ...
                'allowedShapes', "circle", 'radiusRange', [5 5], ...
                'numTargetsRange', [2 2]), 'randomSeed', 42);
            rngBefore = rng;
            randomA = build_cfg(setup);
            testCase.verifyEqual(rng, rngBefore);
            randomB = build_cfg(setup);
            testCase.verifyEqual(randomA.targets, randomB.targets);
            testCase.verifyFalse(any( ...
                randomA.targets(1).mask(:) & randomA.targets(2).mask(:)));

            defaultRandom = build_cfg(struct( ...
                'targetOpts', struct(), 'randomSeed', 11));
            testCase.verifyEqual(numel(defaultRandom.targets), 1);
            testCase.verifyEqual(defaultRandom.targets(1).name, "rectangle");
        end

        function exactTargetsAndOverlapRules(testCase)
            restoreState = configureHybridBuilder(); %#ok<NASGU>
            exact = struct('name', "circle", ...
                'properties', struct('center', [170 200], 'radius', 10), ...
                'material', struct('epsr', 3.7, 'cond_e', 0.04));
            fixed = build_cfg(struct('targetSpecs', exact));
            testCase.verifyEqual(numel(fixed.targets), 1);
            testCase.verifyEqual(fixed.targets(1).material.murx, 1);
            testCase.verifyEqual(fixed.targets(1).material.mury, 1);
            testCase.verifyEqual(fixed.targets(1).material.cond_m, 0);

            opts = struct('allowedShapes', "circle", ...
                'radiusRange', [5 5], 'numTargetsRange', [2 2]);
            combined = build_cfg(struct( ...
                'targetSpecs', exact, 'targetOpts', opts, 'randomSeed', 42));
            testCase.verifyEqual(numel(combined.targets), 3);
            testCase.verifyEqual(combined.targets(1).mask, ...
                fixed.targets(1).mask);
            for targetIdx = 2:3
                testCase.verifyFalse(any(combined.targets(1).mask(:) & ...
                    combined.targets(targetIdx).mask(:)));
            end

            testCase.verifyError(@() build_cfg(struct( ...
                'targetSpecs', [exact exact])), ...
                'build_cfg:OverlappingExactTargets');
            outside = exact;
            outside.properties.center = [300 300];
            testCase.verifyError(@() build_cfg(struct( ...
                'targetSpecs', outside)), 'build_cfg:TargetOutsideDOI');

            large = exact;
            large.properties.center = [200 200];
            large.properties.radius = 70;
            overlapOpts = struct('allowedShapes', "rectangle", ...
                'sideRange', [31 31], 'allowOverlap', true);
            overlapping = build_cfg(struct( ...
                'targetSpecs', large, 'targetOpts', overlapOpts, ...
                'randomSeed', 7));
            overlap = overlapping.targets(1).mask & ...
                overlapping.targets(2).mask;
            testCase.verifyTrue(any(overlap(:)));
            testCase.verifyEqual(overlapping.grid.epsr(overlap), ...
                repmat(overlapping.targets(2).material.epsr, nnz(overlap), 1));
        end

        function inhomogeneityIsMaskedAndReportsOffsets(testCase)
            restoreState = configureHybridBuilder(); %#ok<NASGU>
            cfg = build_cfg();
            mask = cfg.targets(1).mask;
            opts = struct('variance', 0.04, ...
                'correlationCells', 6, 'seed', 123);
            rngBefore = rng;
            [variedA, infoA] = fdtdmat.applyEpsrInhomogeneity( ...
                cfg, mask, opts);
            testCase.verifyEqual(rng, rngBefore);
            [variedB, infoB] = fdtdmat.applyEpsrInhomogeneity( ...
                cfg, mask, opts);
            testCase.verifyEqual(variedA.grid.epsr, variedB.grid.epsr);
            testCase.verifyEqual(infoA, infoB);
            testCase.verifyEqual(variedA.grid.epsr(~mask), ...
                cfg.grid.epsr(~mask));
            offset = variedA.grid.epsr(mask) - cfg.grid.epsr(mask);
            testCase.verifyEqual(infoA.meanOffset, mean(offset), ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(infoA.realizedVarianceOfOffset, ...
                var(offset), 'AbsTol', 1e-12);
            testCase.verifyEqual(infoA.modifiedCellCount, nnz(offset ~= 0));
            testCase.verifyEqual(infoA.minEpsr, ...
                min(variedA.grid.epsr(mask)));
            testCase.verifyEqual(infoA.maxEpsr, ...
                max(variedA.grid.epsr(mask)));
            testCase.verifyEqual(cfg.grid.epsr(mask), ...
                repmat(5, nnz(mask), 1));
        end
    end
end

function cleanup = configureHybridBuilder()
previousPath = path;
previousRng = rng;
cleanup = onCleanup(@() restoreState(previousPath, previousRng));
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'inverse_problem', 'algorithms', ...
    'hybrid_tr_mwi'), '-begin');
end

function restoreState(previousPath, previousRng)
path(previousPath);
rng(previousRng);
end

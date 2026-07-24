classdef testGeometryLibrary < matlab.unittest.TestCase
    %testGeometryLibrary Unit tests for the MATLAB FDTD geometry library.

    methods (TestClassSetup)
        function addLibraryToPath(testCase)
            sourceFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(sourceFolder));
        end
    end

    methods (Test)
        function createGridInitializesCoordinatesAndBackground(testCase)
            background = struct("epsr", 2.5, "murx", 1.2, "mury", 1.3, "cond_e", 0.03, "cond_m", 0.04);

            grid = fdtdmat.createGrid([4 3], [0.1 0.2], background, [1.0 2.0]);

            testCase.verifyEqual(grid.sizeXY, [4 3]);
            testCase.verifyEqual(grid.spacingXY, [0.1 0.2], AbsTol=1e-12);
            testCase.verifyEqual(grid.originPhysical, [1.0 2.0], AbsTol=1e-12);
            testCase.verifyEqual(grid.xIndex(4, 3), 4);
            testCase.verifyEqual(grid.yIndex(4, 3), 3);
            testCase.verifyEqual(grid.xPhysical(4, 3), 1.3, AbsTol=1e-12);
            testCase.verifyEqual(grid.yPhysical(4, 3), 2.4, AbsTol=1e-12);
            testCase.verifyEqual(grid.epsr, 2.5 * ones(4, 3), AbsTol=1e-12);
            testCase.verifyEqual(size(grid.murx), [4 2]);
            testCase.verifyEqual(size(grid.mury), [3 3]);
            testCase.verifyEqual(grid.murx, 1.2 * ones(4, 2), AbsTol=1e-12);
            testCase.verifyEqual(grid.mury, 1.3 * ones(3, 3), AbsTol=1e-12);
            testCase.verifyEqual(grid.cond_e, 0.03 * ones(4, 3), AbsTol=1e-12);
            testCase.verifyEqual(grid.cond_m, 0.04 * ones(4, 3), AbsTol=1e-12);
            testCase.verifyEqual(grid.epsr_bg, grid.epsr, AbsTol=1e-12);
            testCase.verifyEqual(grid.murx_bg, grid.murx, AbsTol=1e-12);
            testCase.verifyEqual(grid.mury_bg, grid.mury, AbsTol=1e-12);
            testCase.verifyEqual(grid.cond_e_bg, grid.cond_e, AbsTol=1e-12);
            testCase.verifyEqual(grid.cond_m_bg, grid.cond_m, AbsTol=1e-12);
        end

        function circleReturnsMaskAndSignedDistance(testCase)
            grid = fdtdmat.createGrid([301 301], [1 1]);

            region = fdtdgeom.shape_circle(grid, [151 151], 18, "index");

            testCase.verifyTrue(region.mask(151, 151));
            testCase.verifyFalse(region.mask(151, 170));
            testCase.verifyEqual(region.signedDistance(151, 151), -18, AbsTol=1e-12);
            testCase.verifyLessThan(region.signedDistance(151, 151), 0);
        end

        function rectangleCreatesExpectedMask(testCase)
            grid = fdtdmat.createGrid([10 10], [1 1]);

            region = fdtdgeom.shape_rectangle(grid, [3 6 4 8], "index");

            testCase.verifyTrue(region.mask(3, 4));
            testCase.verifyTrue(region.mask(5, 6));
            testCase.verifyFalse(region.mask(2, 4));
            testCase.verifyFalse(region.mask(7, 8));
        end

        function booleanRegionOperationsCreateExpectedMasks(testCase)
            grid = fdtdmat.createGrid([5 5], [1 1]);
            left = fdtdgeom.shape_rectangle(grid, [1 3 1 5], "index");
            right = fdtdgeom.shape_rectangle(grid, [3 5 1 5], "index");
            outer = fdtdgeom.shape_rectangle(grid, [1 5 1 5], "index");
            cutout = fdtdgeom.shape_rectangle(grid, [2 4 2 4], "index");
            expectedUnion = true(5, 5);
            expectedIntersection = false(5, 5);
            expectedIntersection(3, :) = true;
            expectedSubtraction = true(5, 5);
            expectedSubtraction(2:4, 2:4) = false;

            unionRegion = fdtdgeom.merge_unionRegions(left, right);
            intersectionRegion = fdtdgeom.merge_intersectRegions(left, right);
            subtractionRegion = fdtdgeom.merge_subtractRegion(outer, cutout);

            testCase.verifyEqual(unionRegion.mask, expectedUnion);
            testCase.verifyEqual(intersectionRegion.mask, expectedIntersection);
            testCase.verifyEqual(subtractionRegion.mask, expectedSubtraction);
        end

        function applyRegionUpdatesOnlyMaskedMaterialCells(testCase)
            grid = fdtdmat.createGrid([20 20], [1 1]);
            region = fdtdgeom.shape_circle(grid, [10 10], 3, "index");

            grid = fdtdmat.applyRegion(grid, region, struct("epsr", 4.0, "murx", 2.0, "mury", 3.0, "cond_e", 0.01));

            testCase.verifyEqual(grid.epsr(10, 10), 4.0, AbsTol=1e-12);
            testCase.verifyEqual(grid.cond_e(10, 10), 0.01, AbsTol=1e-12);
            testCase.verifyEqual(grid.murx(10, 10), 2.0, AbsTol=1e-12);
            testCase.verifyEqual(grid.mury(10, 10), 3.0, AbsTol=1e-12);
            testCase.verifyEqual(grid.epsr(1, 1), 1.0, AbsTol=1e-12);
            testCase.verifyEqual(grid.cond_e(1, 1), 0.0, AbsTol=1e-12);
            testCase.verifyEqual(grid.murx(1, 1), 1.0, AbsTol=1e-12);
            testCase.verifyEqual(grid.mury(1, 1), 1.0, AbsTol=1e-12);
            testCase.verifyEqual(grid.epsr_bg, ones(20, 20), AbsTol=1e-12);
            testCase.verifyEqual(grid.cond_e_bg, zeros(20, 20), AbsTol=1e-12);
        end

        function setBackgroundDefaultPreservesMexContract(testCase)
            cfg = struct('Nx', 4, 'Ny', 3, 'dx', 0.1, 'dy', 0.2, 'marker', 17);
            cfg.grid = struct('originPhysical', [1.0 2.0]);
            background = struct('epsr', 2.5, 'murx', 1.2, 'mury', 1.3, ...
                'cond_e', 0.03, 'cond_m', 0.04);

            cfg = fdtdmat.setBackgroundDefault(cfg, background);

            testCase.verifyEqual(cfg.marker, 17);
            testCase.verifyEqual(cfg.grid.sizeXY, [4 3]);
            testCase.verifyEqual(cfg.grid.spacingXY, [0.1 0.2], AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.originPhysical, [1.0 2.0], AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.epsr, 2.5 * ones(4, 3), AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.murx, 1.2 * ones(4, 2), AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.mury, 1.3 * ones(3, 3), AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.cond_e, 0.03 * ones(4, 3), AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.cond_m, 0.04 * ones(4, 3), AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.epsr_bg, cfg.grid.epsr, AbsTol=1e-12);
            testCase.verifyEqual(cfg.grid.background, background);
        end

        function rectangularPmlMatchesInlineFormula(testCase)
            grid = fdtdmat.createGrid([8 6], [0.1 0.1]);
            thickness = 2;
            order = 3;
            sigmaMax = 10;
            kappaMax = 4;
            [x0, y0] = ndgrid(0:7, 0:5);
            expectedDepthX = max(0, thickness - min(x0, 7 - x0));
            expectedDepthY = max(0, thickness - min(y0, 5 - y0));
            expectedProfileX = (expectedDepthX / thickness) .^ order;
            expectedProfileY = (expectedDepthY / thickness) .^ order;

            pml = fdtdpml.build_rectangularPML(grid, thickness, order, sigmaMax, kappaMax);

            testCase.verifyEqual(pml.depthX, expectedDepthX, AbsTol=1e-12);
            testCase.verifyEqual(pml.depthY, expectedDepthY, AbsTol=1e-12);
            testCase.verifyEqual(pml.profileX, expectedProfileX, AbsTol=1e-12);
            testCase.verifyEqual(pml.profileY, expectedProfileY, AbsTol=1e-12);
            testCase.verifyEqual(pml.condx, sigmaMax * expectedProfileX, AbsTol=1e-12);
            testCase.verifyEqual(pml.condy, sigmaMax * expectedProfileY, AbsTol=1e-12);
            testCase.verifyEqual(pml.kx, 1 + (kappaMax - 1) * expectedProfileX, AbsTol=1e-12);
            testCase.verifyEqual(pml.ky, 1 + (kappaMax - 1) * expectedProfileY, AbsTol=1e-12);
        end

        function radialPmlProjectsRadialProfileOntoCartesianMaps(testCase)
            grid = fdtdmat.createGrid([5 5], [1 1]);
            thickness = 1;
            order = 1;
            sigmaMax = 8;
            kappaMax = 5;
            diagonalProfile = sqrt(2) - 1;

            pml = fdtdpml.build_radialPML(grid, [3 3], thickness, order, sigmaMax, kappaMax, "index");

            testCase.verifyEqual(pml.outerRadius, 2, AbsTol=1e-12);
            testCase.verifyEqual(pml.innerRadius, 1, AbsTol=1e-12);
            testCase.verifyTrue(pml.backgroundMask(3, 3));
            testCase.verifyTrue(pml.innerMask(3, 3));
            testCase.verifyTrue(pml.mask(5, 3));
            testCase.verifyFalse(pml.mask(5, 5));
            testCase.verifyEqual(pml.depth(5, 5), 0, AbsTol=1e-12);
            testCase.verifyEqual(pml.depth(5, 3), 1, AbsTol=1e-12);
            testCase.verifyEqual(pml.condx(5, 3), 8, AbsTol=1e-12);
            testCase.verifyEqual(pml.condy(5, 3), 0, AbsTol=1e-12);
            testCase.verifyEqual(pml.kx(5, 3), 5, AbsTol=1e-12);
            testCase.verifyEqual(pml.ky(5, 3), 1, AbsTol=1e-12);
            testCase.verifyEqual(pml.condx(3, 5), 0, AbsTol=1e-12);
            testCase.verifyEqual(pml.condy(3, 5), 8, AbsTol=1e-12);
            testCase.verifyEqual(pml.condx(4, 4), sigmaMax * diagonalProfile * 0.5, AbsTol=1e-12);
            testCase.verifyEqual(pml.condy(4, 4), sigmaMax * diagonalProfile * 0.5, AbsTol=1e-12);
            testCase.verifyEqual(pml.condx(3, 3), 0, AbsTol=1e-12);
            testCase.verifyEqual(pml.condy(3, 3), 0, AbsTol=1e-12);
        end

        function radialPmlInfersOuterCircleFromNearestGridEdge(testCase)
            grid = fdtdmat.createGrid([251 269], [1 1]);
            thickness = 25;
            order = 1;
            sigmaMax = 8;
            kappaMax = 5;

            pml = fdtdpml.build_radialPML(grid, [126 135], thickness, order, sigmaMax, kappaMax, "index");

            testCase.verifyEqual(pml.outerRadius, 125, AbsTol=1e-12);
            testCase.verifyEqual(pml.innerRadius, 100, AbsTol=1e-12);
            testCase.verifyEqual(pml.depth(226, 135), 0, AbsTol=1e-12);
            testCase.verifyEqual(pml.depth(251, 135), 25, AbsTol=1e-12);
            testCase.verifyEqual(pml.depth(126, 260), 25, AbsTol=1e-12);
            testCase.verifyEqual(pml.profile(241, 135), 0.6, AbsTol=1e-12);
            testCase.verifyEqual(pml.profile(195, 227), 0.6, AbsTol=1e-12);
            testCase.verifyFalse(pml.backgroundMask(126, 269));
            testCase.verifyFalse(pml.mask(126, 269));
            testCase.verifyEqual(pml.depth(126, 269), 0, AbsTol=1e-12);
            testCase.verifyEqual(pml.condr(126, 269), 0, AbsTol=1e-12);
            testCase.verifyEqual(pml.kappar(126, 269), 1, AbsTol=1e-12);
        end
        function driverPreservesCsvContract(testCase)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            exampleDir = fullfile(repoRoot, 'forward_solver', 'examples');

            scriptText = [fileread(fullfile(exampleDir, 'build_fdtd_case.m')), newline, ...
                fileread(fullfile(exampleDir, 'run_fdtd_file_case.m'))];

            testCase.verifyNotEmpty(strfind(scriptText, 'epsr.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'murx.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'mury.csv'));
            legacyName = ['m' 'u' 'r' '.' 'csv'];
            testCase.verifyEmpty(strfind(scriptText, legacyName));
            testCase.verifyNotEmpty(strfind(scriptText, 'cond_e.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'cond_m.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'condx.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'condy.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'kx.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'ky.csv'));
            testCase.verifyNotEmpty(strfind(scriptText, 'fdtdmat.createGrid'));
            testCase.verifyNotEmpty(strfind(scriptText, 'fdtdpml.build_radialPML'));
            testCase.verifyEmpty(strfind(scriptText, 'pmlRadius'));
            testCase.verifyNotEmpty(strfind(scriptText, 'etaBackground'));
            testCase.verifyEmpty(regexp(scriptText, '\n\s*maps\s*=', 'once'));
            testCase.verifyEmpty(regexp(scriptText, '\n\s*maps\.', 'once'));
        end
    end
end

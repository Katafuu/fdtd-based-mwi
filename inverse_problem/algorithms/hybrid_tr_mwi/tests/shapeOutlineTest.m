classdef shapeOutlineTest < matlab.unittest.TestCase
    %shapeOutlineTest Tests shape-point filtering, ordering, and plotting.

    methods (TestClassSetup)
        function addHybridPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            hybridDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(hybridDir, ...
                IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testCleanPointCloudRetainsEveryPoint(testCase)
            points = [0 0; 1 0; 1 1; 0 1];

            [inliers, outliers, inlierIndices, outlierIndices] = ...
                filterShapeOutlinePoints(points, 2);

            testCase.verifyEqual(inliers, points);
            testCase.verifyEmpty(outliers);
            testCase.verifyEqual(inlierIndices, (1:4).');
            testCase.verifyEmpty(outlierIndices);
        end

        function testExtremeXAndYCoordinatesAreExcluded(testCase)
            points = [ ...
                0 0; 1 0; 0 1; 1 1; -1 0; 0 -1; ...
                -1 -1; 1 -1; -1 1; 0.5 0.5; 100 0; 0 100];

            [~, outliers, ~, outlierIndices] = ...
                filterShapeOutlinePoints(points, 2);

            testCase.verifyEqual(outliers, [100 0; 0 100]);
            testCase.verifyEqual(outlierIndices, [11; 12]);
        end

        function testNonfiniteCoordinatesAreExcluded(testCase)
            points = [0 0; NaN 1; 1 Inf; 1 1];

            [inliers, outliers, inlierIndices, outlierIndices] = ...
                filterShapeOutlinePoints(points, 2);

            testCase.verifyEqual(inliers, [0 0; 1 1]);
            testCase.verifyEqual(outliers, [NaN 1; 1 Inf]);
            testCase.verifyEqual(inlierIndices, [1; 4]);
            testCase.verifyEqual(outlierIndices, [2; 3]);
        end

        function testConstantCoordinateDoesNotCreateOutliers(testCase)
            points = [5 0; 5 1; 5 2; 5 3; 5 4];

            [inliers, outliers, ~, ~, statistics] = ...
                filterShapeOutlinePoints(points, 2);

            testCase.verifyEqual(inliers, points);
            testCase.verifyEmpty(outliers);
            testCase.verifyEqual(statistics.StandardDeviation(1), 0, ...
                AbsTol=eps);
        end

        function testDoiFilteringPrecedesStatisticalFiltering(testCase)
            doiMask = false(20, 20);
            doiMask(5:15, 5:15) = true;
            points = [ ...
                9 9; 10 9; 11 9; 9 10; 10 10; ...
                11 10; 9 11; 10 11; 14 14; 100 100];

            [doiPoints, doiIndices, outsidePoints, outsideIndices] = ...
                selectShapePointsInsideDoi(points, doiMask);
            [~, statisticalOutliers] = filterShapeOutlinePoints(doiPoints);

            testCase.verifyEqual(doiPoints, points(1:9, :));
            testCase.verifyEqual(doiIndices, (1:9).');
            testCase.verifyEqual(outsidePoints, [100 100]);
            testCase.verifyEqual(outsideIndices, 10);
            testCase.verifyEqual(statisticalOutliers, [14 14]);
        end

        function testNearestTourUsesDeterministicTieBreak(testCase)
            points = [0 0; 1 0; 0 1; 1 1];
            originalIndices = [10; 20; 30; 40];

            [orderedPoints, orderedIndices] = ...
                orderShapeOutlinePoints(points, originalIndices);

            testCase.verifyEqual(orderedIndices, [10; 20; 40; 30]);
            testCase.verifyEqual(orderedPoints, points([1 2 4 3], :));
        end

        function testClosedTourUsesEveryPointWithDegreeTwo(testCase)
            points = [0 0; 2 0; 2 2; 0 2; 1 1];
            originalIndices = [2; 4; 6; 8; 10];

            [~, orderedIndices] = ...
                orderShapeOutlinePoints(points, originalIndices);
            closedIndices = [orderedIndices; orderedIndices(1)];
            edgeVertices = [closedIndices(1:end-1); closedIndices(2:end)];

            testCase.verifyEqual(sort(orderedIndices), sort(originalIndices));
            testCase.verifyEqual(sort(edgeVertices), ...
                sort(repelem(originalIndices, 2)));
            testCase.verifyEqual(closedIndices(1), closedIndices(end));
        end

        function testSplitOverviewsSeparatePointsAndOutline(testCase)
            [cfg, targetMask] = shapeOutlineTest.createPlotInputs();
            points = [ ...
                9 9; 10 9; 11 9; 11 10; 11 11; ...
                10 11; 9 11; 9 10; 10 10; 20 20; 29 29];

            [pointsFigure, outlineFigure] = plotShapeEstimateOverview( ...
                points, cfg, targetMask, 1.5);
            testCase.addTeardown(@() close([pointsFigure outlineFigure]));
            pointsTarget = findobj(pointsFigure, 'Tag', 'targetOutline');
            outlineTarget = findobj(outlineFigure, 'Tag', 'targetOutline');
            inlierHandle = findobj(pointsFigure, 'Tag', 'inlierPoints');
            outlierHandle = findobj(pointsFigure, 'Tag', 'outlierPoints');
            pointsEstimate = findobj(pointsFigure, 'Tag', 'estimatedOutline');
            outlineEstimate = findobj(outlineFigure, 'Tag', 'estimatedOutline');
            outlineInliers = findobj(outlineFigure, 'Tag', 'inlierPoints');
            outlineOutliers = findobj(outlineFigure, 'Tag', 'outlierPoints');

            testCase.verifyNumElements(pointsTarget, 1);
            testCase.verifyNumElements(outlineTarget, 1);
            testCase.verifyNumElements(inlierHandle, 1);
            testCase.verifyNumElements(outlierHandle, 1);
            testCase.verifyEmpty(pointsEstimate);
            testCase.verifyNumElements(outlineEstimate, 1);
            testCase.verifyEmpty(outlineInliers);
            testCase.verifyEmpty(outlineOutliers);
            testCase.verifyEqual( ...
                [outlineEstimate.XData(1) outlineEstimate.YData(1)], ...
                [outlineEstimate.XData(end) outlineEstimate.YData(end)]);
            testCase.verifyEqual( ...
                [outlierHandle.XData(:) outlierHandle.YData(:)], [20 20]);
        end

        function testFewerThanThreeInliersWarns(testCase)
            [cfg, targetMask] = shapeOutlineTest.createPlotInputs();
            points = [10 10; 12 12];
            testCase.addTeardown(@shapeOutlineTest.closeOverviewFigures);

            testCase.verifyWarning( ...
                @() plotShapeEstimateOverview(points, cfg, targetMask, 1.5), ...
                'plotShapeEstimateOverview:InsufficientInliers');
        end
    end

    methods (Static, Access = private)
        function [cfg, targetMask] = createPlotInputs()
            cfg = struct();
            cfg.Nx = 30;
            cfg.Ny = 30;
            cfg.antennas.doiMask = false(cfg.Nx, cfg.Ny);
            cfg.antennas.doiMask(4:27, 4:27) = true;
            cfg.antennas.pos = [4 4; 27 4; 27 27; 4 27];
            targetMask = false(cfg.Nx, cfg.Ny);
            targetMask(12:18, 10:20) = true;
        end

        function closeOverviewFigures()
            close(findall(groot, 'Type', 'figure', ...
                '-regexp', 'Tag', '^shapeEstimate'));
        end
    end
end

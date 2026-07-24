classdef plotCircularFootprintCoverageTest < matlab.unittest.TestCase
    %plotCircularFootprintCoverageTest Tests the hybrid mask diagnostic.

    properties
        Axes
        Coverage
        ExpectedCoverage
        FocusPoint
        AntennaPositions
        XRange
        YRange
        Cfg
    end

    methods (TestClassSetup)
        function addHybridPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            hybridDir = fileparts(fileparts(mfilename('fullpath')));
            algorithmsDir = fileparts(hybridDir);
            testCase.applyFixture(PathFixture(hybridDir, ...
                IncludingSubfolders=true));
            testCase.applyFixture(PathFixture(fullfile(algorithmsDir, ...
                'time_reversal', 'lib')));
        end
    end

    methods (TestMethodSetup)
        function createCoveragePlot(testCase)
            xRange = 10:14;
            yRange = 20:25;
            mask4D = false(numel(xRange), numel(yRange), 3, 3);
            mask4D(1, 1, 1, 1) = true;
            mask4D(2:3, 2:3, 1, 2) = true;
            mask4D(3:4, 3:4, 2, 1) = true;
            mask4D(5, 6, 1, 3) = true;

            pairAccepted = false(3, 3);
            pairAccepted(1, 2) = true;
            pairAccepted(2, 1) = true;

            cfg = struct();
            cfg.Nx = 30;
            cfg.Ny = 35;
            cfg.antennas.pos = [10 20; 12 22; 14 25];
            cfg.antennas.txAntennas = 2;
            cfg.antennas.doiMask = false(cfg.Nx, cfg.Ny);
            cfg.antennas.doiMask(8:16, 18:27) = true;
            cfg.targets = struct('name', 'circle', 'properties', ...
                struct('center', [12 23], 'radius', 2));
            focusPoint = [12 23];

            expectedCoverage = zeros(cfg.Nx, cfg.Ny);
            expectedCoverage(xRange, yRange) = ...
                double(mask4D(:, :, 1, 2)) + double(mask4D(:, :, 2, 1));
            figureHandle = figure('Visible', 'off');
            testCase.addTeardown(@() close(figureHandle));
            axesHandle = axes('Parent', figureHandle);
            axes(axesHandle);
            coverage = plotCircularFootprintCoverage( ...
                mask4D, xRange, yRange, pairAccepted, cfg, focusPoint);

            testCase.Axes = axesHandle;
            testCase.Coverage = coverage;
            testCase.ExpectedCoverage = expectedCoverage;
            testCase.FocusPoint = focusPoint;
            testCase.AntennaPositions = cfg.antennas.pos;
            testCase.XRange = xRange;
            testCase.YRange = yRange;
            testCase.Cfg = cfg;
        end
    end

    methods (Test)
        function testCoverageCountsAcceptedDirectedMasks(testCase)
            testCase.verifyEqual(testCase.Coverage, ...
                testCase.ExpectedCoverage);
        end

        function testReturnedAndDisplayedCoverageUseFullGrid(testCase)
            imageHandle = findobj(testCase.Axes, 'Tag', 'maskCoverage');
            testCase.verifySize(testCase.Coverage, [30 35]);
            testCase.verifySize(imageHandle.CData, [35 30]);
        end

        function testRejectedAndSelfMasksAreExcluded(testCase)
            testCase.verifyEqual(testCase.Coverage(10, 20), 0);
            testCase.verifyEqual(testCase.Coverage(14, 25), 0);
        end

        function testCoverageIsZeroOutsideInteriorRanges(testCase)
            outsideCoverage = testCase.Coverage;
            outsideCoverage(testCase.XRange, testCase.YRange) = 0;
            testCase.verifyEqual(outsideCoverage, ...
                zeros(testCase.Cfg.Nx, testCase.Cfg.Ny));
        end

        function testUsesExactFullGridAxesLimits(testCase)
            testCase.verifyEqual(testCase.Axes.XLim, [0.5 30.5]);
            testCase.verifyEqual(testCase.Axes.YLim, [0.5 35.5]);
            testCase.verifyEqual(testCase.Axes.YDir, 'normal');
            testCase.verifyEqual(testCase.Axes.DataAspectRatio, [1 1 1]);
        end

        function testReciprocalPairHasSingleOutline(testCase)
            outlineHandles = findall(testCase.Axes, 'Tag', 'maskOutline');
            testCase.verifyEqual(numel(outlineHandles), 1);
        end

        function testCoverageUsesXYDisplayOrientation(testCase)
            imageHandle = findobj(testCase.Axes, 'Tag', 'maskCoverage');
            testCase.verifyEqual(numel(imageHandle), 1);
            testCase.verifyEqual(imageHandle.CData, ...
                testCase.ExpectedCoverage.');
        end

        function testFocusAndAntennaOverlaysUseGridCoordinates(testCase)
            focusHandle = findobj(testCase.Axes, 'Tag', 'trFocus');
            antennaHandle = findobj(testCase.Axes, 'Tag', 'antennaPositions');
            testCase.verifyEqual([focusHandle.XData focusHandle.YData], ...
                testCase.FocusPoint);
            testCase.verifyEqual([antennaHandle.XData(:) antennaHandle.YData(:)], ...
                testCase.AntennaPositions);
        end

        function testAllDomainAndMaskGraphicsExist(testCase)
            targetHandles = findall(testCase.Axes, 'Tag', 'targetOutline');
            antennaHandles = findall(testCase.Axes, 'Tag', 'antennaPositions');
            txHandles = findall(testCase.Axes, 'Tag', 'txAntenna');
            doiHandles = findall(testCase.Axes, 'Tag', 'doiOutline');
            maskHandles = findall(testCase.Axes, 'Tag', 'maskOutline');
            focusHandles = findall(testCase.Axes, 'Tag', 'trFocus');

            testCase.verifyNumElements(targetHandles, 1);
            testCase.verifyNumElements(antennaHandles, 1);
            testCase.verifyNumElements(txHandles, 1);
            testCase.verifyNumElements(doiHandles, 1);
            testCase.verifyNumElements(maskHandles, 1);
            testCase.verifyNumElements(focusHandles, 1);
        end
    end
end
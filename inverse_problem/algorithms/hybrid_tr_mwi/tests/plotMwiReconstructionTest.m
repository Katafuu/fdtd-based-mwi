classdef plotMwiReconstructionTest < matlab.unittest.TestCase
    %plotMwiReconstructionTest Tests the full-grid reconstruction renderer.

    properties
        Axes
        Cfg
        EpsrAvgFull
        Image
        XRange
        YRange
    end

    methods (TestClassSetup)
        function addPlottingPaths(testCase)
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
        function createReconstructionPlot(testCase)
            cfg = struct();
            cfg.Nx = 18;
            cfg.Ny = 14;
            cfg.grid.background.epsr = 1.7;
            cfg.antennas.pos = [3 7; 9 2; 16 7; 9 13];
            cfg.antennas.txAntennas = 2;
            cfg.antennas.doiMask = false(cfg.Nx, cfg.Ny);
            cfg.antennas.doiMask(4:15, 3:12) = true;

            propertiesTemplate = struct('center', [], 'radius', [], ...
                'bounds', [], 'vertices', []);
            circleProperties = propertiesTemplate;
            circleProperties.center = [7 7];
            circleProperties.radius = 2;
            rectangleProperties = propertiesTemplate;
            rectangleProperties.bounds = [10 13 4 7];
            triangleProperties = propertiesTemplate;
            triangleProperties.vertices = [8 9; 11 11; 6 12];
            cfg.targets = [ ...
                struct('name', 'circle', 'properties', circleProperties), ...
                struct('name', 'rectangle', 'properties', rectangleProperties), ...
                struct('name', 'triangle', 'properties', triangleProperties)];

            xRange = 4:15;
            yRange = 3:12;
            epsrAvgFull = cfg.grid.background.epsr .* ...
                ones(cfg.Nx, cfg.Ny);
            epsrAvgFull(xRange, yRange) = reshape( ...
                linspace(2, 4, numel(xRange) * numel(yRange)), ...
                numel(xRange), numel(yRange));

            figureHandle = figure('Visible', 'off');
            testCase.addTeardown(@() close(figureHandle));
            axesHandle = axes('Parent', figureHandle);
            imageHandle = plotMwiReconstruction( ...
                axesHandle, cfg, epsrAvgFull);

            testCase.Axes = axesHandle;
            testCase.Cfg = cfg;
            testCase.EpsrAvgFull = epsrAvgFull;
            testCase.Image = imageHandle;
            testCase.XRange = xRange;
            testCase.YRange = yRange;
        end
    end

    methods (Test)
        function testFullGridOrientationAndLimits(testCase)
            testCase.verifySize(testCase.Image.CData, ...
                [testCase.Cfg.Ny testCase.Cfg.Nx]);
            testCase.verifyEqual(testCase.Image.CData, ...
                testCase.EpsrAvgFull.');
            testCase.verifyEqual(testCase.Axes.XLim, [0.5 18.5]);
            testCase.verifyEqual(testCase.Axes.YLim, [0.5 14.5]);
            testCase.verifyEqual(testCase.Axes.YDir, 'normal');
            testCase.verifyEqual(testCase.Axes.DataAspectRatio, [1 1 1]);
        end

        function testDisplaysConfiguredBackgroundPadding(testCase)
            displayedData = testCase.Image.CData.';
            displayedData(testCase.XRange, testCase.YRange) = ...
                testCase.Cfg.grid.background.epsr;
            testCase.verifyEqual(displayedData, ...
                testCase.Cfg.grid.background.epsr .* ...
                ones(testCase.Cfg.Nx, testCase.Cfg.Ny));
        end

        function testDrawsAllSupportedTargetOutlines(testCase)
            targetHandles = findall(testCase.Axes, 'Tag', 'targetOutline');
            pointCounts = sort(arrayfun( ...
                @(handle) numel(handle.XData), targetHandles));
            testCase.verifyNumElements(targetHandles, 3);
            testCase.verifyEqual(pointCounts, [4; 5; 256]);
        end

        function testDrawsAntennaTransmitterAndDoi(testCase)
            antennaHandle = findall(testCase.Axes, 'Tag', 'antennaPositions');
            txHandle = findall(testCase.Axes, 'Tag', 'txAntenna');
            doiHandle = findall(testCase.Axes, 'Tag', 'doiOutline');
            testCase.verifyNumElements(antennaHandle, 1);
            testCase.verifyNumElements(txHandle, 1);
            testCase.verifyNumElements(doiHandle, 1);
            testCase.verifyEqual([antennaHandle.XData(:) ...
                antennaHandle.YData(:)], testCase.Cfg.antennas.pos);
            testCase.verifyEqual([txHandle.XData txHandle.YData], ...
                testCase.Cfg.antennas.pos(2, :));
            testCase.verifyEqual(antennaHandle.Color, [1 1 1]);
        end
    end
end

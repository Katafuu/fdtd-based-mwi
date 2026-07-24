classdef plot_mwiTest < matlab.unittest.TestCase
    %plot_mwiTest Tests six-figure hybrid plotting orchestration.

    methods (TestClassSetup)
        function addPlottingPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            hybridDir = fileparts(fileparts(mfilename('fullpath')));
            algorithmsDir = fileparts(hybridDir);
            testCase.applyFixture(PathFixture(hybridDir, ...
                IncludingSubfolders=true));
            testCase.applyFixture(PathFixture(fullfile(algorithmsDir, ...
                'time_reversal', 'lib')));
            testCase.applyFixture(PathFixture(fullfile(algorithmsDir, ...
                'transmission_mwi', 'lib')));
        end
    end

    methods (Test)
        function testExportsFilesDimensionsDpiAndRestoresState(testCase)
            [cfg, plotData] = plot_mwiTest.createPlotInputs();
            outputDir = string(tempname);
            mkdir(outputDir);
            testCase.addTeardown(@() rmdir(outputDir, 's'));

            originalVisibility = get(groot, 'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot, ...
                'DefaultFigureVisible', originalVisibility));
            set(groot, 'DefaultFigureVisible', 'on');
            sentinelFigure = figure('Visible', 'off');
            testCase.addTeardown(@() close(sentinelFigure));
            figuresBefore = findall(groot, 'Type', 'figure');

            plot_mwi(cfg, plotData, outputDir);

            expectedNames = sort([ ...
                "01_antenna_array.png", "02_tr_actual_target.png", ...
                "03_tr_max_magnitude_focus.png", ...
                "04_tr_minimum_r_focus.png", ...
                "05_footprint_coverage.png", ...
                "06_epsr_reconstruction.png"]);
            pngFiles = dir(fullfile(outputDir, '*.png'));
            actualNames = sort(string({pngFiles.name}));
            pngPaths = cellstr(fullfile(outputDir, actualNames));
            pngInfo = cellfun(@imfinfo, pngPaths, 'UniformOutput', false);
            maxInfo = pngInfo{actualNames == ...
                "03_tr_max_magnitude_focus.png"};
            footprintInfo = pngInfo{actualNames == ...
                "05_footprint_coverage.png"};
            reconstructionInfo = pngInfo{actualNames == ...
                "06_epsr_reconstruction.png"};
            xDpi = cellfun(@(info) info.XResolution, pngInfo) .* 0.0254;
            yDpi = cellfun(@(info) info.YResolution, pngInfo) .* 0.0254;
            figuresAfter = findall(groot, 'Type', 'figure');

            testCase.verifyEqual(actualNames, expectedNames);
            testCase.verifyTrue(all([pngFiles.bytes] > 0));
            testCase.verifyEqual([footprintInfo.Width footprintInfo.Height], ...
                [maxInfo.Width maxInfo.Height]);
            testCase.verifyEqual( ...
                [reconstructionInfo.Width reconstructionInfo.Height], ...
                [maxInfo.Width maxInfo.Height]);
            testCase.verifyEqual(xDpi, 300 .* ones(size(xDpi)), ...
                AbsTol=0.1);
            testCase.verifyEqual(yDpi, 300 .* ones(size(yDpi)), ...
                AbsTol=0.1);
            testCase.verifyEqual(string(get(groot, ...
                'DefaultFigureVisible')), "on");
            testCase.verifyEmpty(setxor( ...
                [figuresBefore.Number], [figuresAfter.Number]));
            testCase.verifyTrue(isgraphics(sentinelFigure, 'figure'));
        end

        function testRestoresStateAfterPlottingFailure(testCase)
            [cfg, plotData] = plot_mwiTest.createPlotInputs();
            plotData.epsrAvgFull = ones(cfg.Nx - 1, cfg.Ny);
            outputDir = string(tempname);
            mkdir(outputDir);
            testCase.addTeardown(@() rmdir(outputDir, 's'));

            originalVisibility = get(groot, 'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot, ...
                'DefaultFigureVisible', originalVisibility));
            set(groot, 'DefaultFigureVisible', 'on');
            sentinelFigure = figure('Visible', 'off');
            testCase.addTeardown(@() close(sentinelFigure));
            figuresBefore = findall(groot, 'Type', 'figure');

            testCase.verifyError(@() plot_mwi(cfg, plotData, outputDir), ...
                'plotMwiReconstruction:GridSizeMismatch');
            figuresAfter = findall(groot, 'Type', 'figure');

            testCase.verifyEqual(string(get(groot, ...
                'DefaultFigureVisible')), "on");
            testCase.verifyEmpty(setxor( ...
                [figuresBefore.Number], [figuresAfter.Number]));
            testCase.verifyTrue(isgraphics(sentinelFigure, 'figure'));
        end
    end

    methods (Static, Access = private)
        function [cfg, plotData] = createPlotInputs()
            cfg = struct();
            cfg.Nx = 24;
            cfg.Ny = 20;
            cfg.pml.thickness = 2;
            cfg.grid.background.epsr = 1.5;
            cfg.grid.epsr = cfg.grid.background.epsr .* ...
                ones(cfg.Nx, cfg.Ny);
            cfg.grid.epsr(10:14, 8:12) = 3;
            cfg.antennas.center = [12 10];
            cfg.antennas.pos = [4 10; 12 3; 21 10; 12 18];
            cfg.antennas.txAntennas = 2;
            cfg.antennas.doiMask = false(cfg.Nx, cfg.Ny);
            cfg.antennas.doiMask(3:22, 3:18) = true;
            cfg.targets = struct('name', 'circle', 'properties', ...
                struct('center', [12 10], 'radius', 3));

            trResult = struct();
            trResult.focus_mag_image = zeros(cfg.Nx, cfg.Ny);
            trResult.focus_mag_image(12, 10) = 4;
            trResult.focus_entropy_image = zeros(cfg.Nx, cfg.Ny);
            trResult.focus_entropy_image(13, 10) = 3;
            trResult.focus_mag_step = 17;
            trResult.focus_entropy_step = 23;

            xRange = 3:22;
            yRange = 3:18;
            mask4D = false(numel(xRange), numel(yRange), 4, 4);
            mask4D(6:14, 5:12, 1, 2) = true;
            mask4D(8:16, 6:13, 2, 1) = true;
            pairAccepted = false(4, 4);
            pairAccepted(1, 2) = true;
            pairAccepted(2, 1) = true;
            epsrAvgFull = cfg.grid.background.epsr .* ...
                ones(cfg.Nx, cfg.Ny);
            epsrAvgFull(xRange, yRange) = 2.5;

            plotData = struct();
            plotData.trResult = trResult;
            plotData.focusPoint = [12 10];
            plotData.mask4D = mask4D;
            plotData.xRange = xRange;
            plotData.yRange = yRange;
            plotData.pairAccepted = pairAccepted;
            plotData.epsrAvgFull = epsrAvgFull;
        end
    end
end

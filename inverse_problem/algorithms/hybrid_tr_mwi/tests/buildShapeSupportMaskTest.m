classdef buildShapeSupportMaskTest < matlab.unittest.TestCase
    %buildShapeSupportMaskTest Verify outline rasterization and DOI clipping.

    methods (TestClassSetup)
        function addHybridLibrary(testCase)
            libraryDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                libraryDir));
        end
    end

    methods (Test)
        function testFilledMaskUsesXYGridOrientation(testCase)
            doiMask = true(18, 14);
            points = [5 4; 11 4; 11 8; 5 8; 100 100];

            [supportMask, outline] = buildShapeSupportMask( ...
                points, doiMask);
            expected = false(18, 14);
            expected(5:11, 4:8) = true;

            testCase.verifyEqual(supportMask, expected);
            testCase.verifySize(outline, [4 2]);
        end

        function testClipsOutlineToDoi(testCase)
            doiMask = true(18, 14);
            doiMask(7:8, 6:7) = false;
            points = [5 4; 11 4; 11 8; 5 8];

            supportMask = buildShapeSupportMask(points, doiMask);
            expected = false(18, 14);
            expected(5:11, 4:8) = true;
            expected = expected & doiMask;

            testCase.verifyEqual(supportMask, expected);
        end

        function testRejectsDegenerateOutline(testCase)
            doiMask = true(18, 14);
            testCase.verifyError( ...
                @() buildShapeSupportMask([5 4; 11 4], doiMask), ...
                'buildShapeSupportMask:InvalidOutline');
            testCase.verifyError( ...
                @() buildShapeSupportMask([5 4; 8 4; 11 4], doiMask), ...
                'buildShapeSupportMask:InvalidOutline');
        end
    end
end

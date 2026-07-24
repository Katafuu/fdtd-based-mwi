classdef buildFootprintMaskTest < matlab.unittest.TestCase
    %buildFootprintMaskTest Tests midpoint-to-focus footprint blending.

    properties (TestParameter)
        invalidBias = struct('belowZero', -0.1, 'aboveOne', 1.1)
    end

    methods (TestClassSetup)
        function addHybridPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            hybridDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(hybridDir, ...
                IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testZeroBiasUsesMidpoint(testCase)
            mask = buildFootprintMaskTest.createMask(0);
            expected = buildFootprintMaskTest.squareMask(4:6, 2:4);

            testCase.verifyEqual(mask(:, :, 1, 1), expected);
        end

        function testUnitBiasUsesFocusProjection(testCase)
            mask = buildFootprintMaskTest.createMask(1);
            expected = buildFootprintMaskTest.squareMask(6:8, 2:4);

            testCase.verifyEqual(mask(:, :, 1, 1), expected);
        end

        function testHalfBiasInterpolatesPosition(testCase)
            mask = buildFootprintMaskTest.createMask(0.5);
            expected = buildFootprintMaskTest.squareMask(5:7, 2:4);

            testCase.verifyEqual(mask(:, :, 1, 1), expected);
        end

        function testOffRayFocusWithinHalfWidthUsesProjection(testCase)
            mask = buildFootprintMaskTest.createConfiguredMask( ...
                1, [7 3.5], true(9, 5));
            expected = buildFootprintMaskTest.squareMask(6:8, 2:4);

            testCase.verifyEqual(mask(:, :, 1, 1), expected);
        end

        function testFocusAtHalfWidthProducesEmptyFootprint(testCase)
            mask = buildFootprintMaskTest.createConfiguredMask( ...
                1, [7 4], true(9, 5));

            testCase.verifyFalse(any(mask, 'all'));
        end

        function testFocusBeyondHalfWidthProducesEmptyFootprint(testCase)
            mask = buildFootprintMaskTest.createConfiguredMask( ...
                1, [7 5], true(9, 5));

            testCase.verifyFalse(any(mask, 'all'));
        end

        function testZeroBiasBypassesDistanceGate(testCase)
            mask = buildFootprintMaskTest.createConfiguredMask( ...
                0, [7 5], true(9, 5));
            expected = buildFootprintMaskTest.squareMask(4:6, 2:4);

            testCase.verifyEqual(mask(:, :, 1, 1), expected);
        end

        function testFootprintIsClippedToDoi(testCase)
            doiMask = true(9, 5);
            doiMask(8, :) = false;
            mask = buildFootprintMaskTest.createConfiguredMask( ...
                1, [7 3], doiMask);
            expected = buildFootprintMaskTest.squareMask(6:7, 2:4);

            testCase.verifyEqual(mask(:, :, 1, 1), expected);
        end

        function testOutOfRangeBiasErrors(testCase, invalidBias)
            testCase.verifyError( ...
                @() buildFootprintMaskTest.createMask(invalidBias), ...
                'buildFootprintMask:InvalidFocusBias');
        end
    end

    methods (Static, Access = private)
        function mask = createMask(focusBias)
            mask = buildFootprintMaskTest.createConfiguredMask( ...
                focusBias, [7 3], true(9, 5));
        end

        function mask = createConfiguredMask(focusBias, focusPoint, doiMask)
            mask = buildFootprintMask(9, 5, 0, [1 3], [9 3], ...
                1, focusPoint, focusBias, doiMask);
        end

        function mask = squareMask(xIndices, yIndices)
            mask = false(9, 5);
            mask(xIndices, yIndices) = true;
        end
    end
end
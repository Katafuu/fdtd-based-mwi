function [footprintMask, xRange, yRange] = ...
        buildFootprintMask(Nx, Ny, pmlThickness, transmitter_cell_indices, receiver_cell_indices, Lfp_cells, focusPoint, focusBias, doiMask)
% buildFootprintMask Construct time-reversal-guided footprint masks.
%
% Antenna index convention:
%   transmitter_cell_indices(:,1) = x index
%   transmitter_cell_indices(:,2) = y index
%   receiver_cell_indices(:,1)    = x index
%   receiver_cell_indices(:,2)    = y index
%
% Matrix convention used here:
%   footprintMask(x, y, tx, rx)
%
% focusPoint is [x y] in the same full-grid index convention. Each
% focusBias blends each centre from the TX-RX midpoint (0) to the point
% on the finite segment nearest focusPoint (1). Focus-guided footprints
% are retained only when that point is less than Lfp_cells from focusPoint.

    validateattributes(focusPoint, {'numeric'}, {'real', 'finite', 'vector', 'numel', 2}, ...
        mfilename, 'focusPoint');
    validateattributes(focusBias, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'focusBias');
    if focusBias < 0 || focusBias > 1
        error('buildFootprintMask:InvalidFocusBias', ...
            'focusBias must be between 0 and 1 inclusive.');
    end
    validateattributes(doiMask, {'logical'}, {'size', [Nx Ny]}, ...
        mfilename, 'doiMask');
    focusPoint = reshape(double(focusPoint), 1, 2);

    Ntx = size(transmitter_cell_indices, 1);
    Nrx = size(receiver_cell_indices, 1);

    xRange = pmlThickness + 1 : Nx - pmlThickness;
    yRange = pmlThickness + 1 : Ny - pmlThickness;

    [Xidx, Yidx] = ndgrid(xRange, yRange);
    croppedDoiMask = doiMask(xRange, yRange);

    footprintMask = false(numel(xRange), numel(yRange), Ntx, Nrx);

    for tx = 1:Ntx
        x1 = transmitter_cell_indices(tx, 1);
        y1 = transmitter_cell_indices(tx, 2);

        for rx = 1:Nrx
            x2 = receiver_cell_indices(rx, 1);
            y2 = receiver_cell_indices(rx, 2);

            % Project the TR focus onto the finite TX-RX segment.
            ray = [x2 - x1, y2 - y1];
            rayLengthSquared = dot(ray, ray);
            if rayLengthSquared == 0
                xip = x1;
                yip = y1;
            else
                tFocus = dot(focusPoint - [x1, y1], ray) / rayLengthSquared;
                tFocus = max(0, min(1, tFocus));
                focusProjection = [x1, y1] + tFocus * ray;
                perpendicularDistance = hypot( ...
                    focusPoint(1) - focusProjection(1), ...
                    focusPoint(2) - focusProjection(2));
                if focusBias > 0 && perpendicularDistance >= Lfp_cells
                    continue;
                end
                t = (1 - focusBias)*0.5 + focusBias*tFocus;
                xip = x1 + t * ray(1);
                yip = y1 + t * ray(2);
            end

            % Square footprint, matching the inequality form in the paper.
            candidateMask = ...
                abs(Xidx - xip) <= Lfp_cells & ...
                abs(Yidx - yip) <= Lfp_cells;
            footprintMask(:, :, tx, rx) = ...
                maskCleanup(candidateMask, croppedDoiMask);
        end
    end
end

function cleanedMask = maskCleanup(candidateMask, croppedDoiMask)
    cleanedMask = candidateMask & croppedDoiMask;
end

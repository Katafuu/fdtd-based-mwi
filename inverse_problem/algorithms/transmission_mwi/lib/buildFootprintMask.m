function [footprintMask, xRange, yRange] = buildFootprintMask(Nx, Ny, pmlThickness, transmitter_cell_indices, receiver_cell_indices, Lfp_cells)
% buildFootprintMask Construct paper-style footprint masks.
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
% The TDS papers map each transmitter-receiver estimate to a localized
% footprint centered at the midpoint/intersection point of the ray, not to
% the entire ray path.

    Ntx = size(transmitter_cell_indices, 1);
    Nrx = size(receiver_cell_indices, 1);

    xRange = pmlThickness + 1 : Nx - pmlThickness;
    yRange = pmlThickness + 1 : Ny - pmlThickness;

    [Xidx, Yidx] = ndgrid(xRange, yRange);

    footprintMask = false(numel(xRange), numel(yRange), Ntx, Nrx);

    for tx = 1:Ntx
        x1 = transmitter_cell_indices(tx, 1);
        y1 = transmitter_cell_indices(tx, 2);

        for rx = 1:Nrx
            x2 = receiver_cell_indices(rx, 1);
            y2 = receiver_cell_indices(rx, 2);

            % Midpoint of the transmitter-receiver path. For the planar
            % simplified geometry, this is the intersection point with the
            % reconstruction plane midway between the antenna plates.
            xip = 0.5 * (x1 + x2);
            yip = 0.5 * (y1 + y2);

            % Square footprint, matching the inequality form in the paper.
            footprintMask(:, :, tx, rx) = ...
                abs(Xidx - xip) <= Lfp_cells & ...
                abs(Yidx - yip) <= Lfp_cells;
        end
    end
end

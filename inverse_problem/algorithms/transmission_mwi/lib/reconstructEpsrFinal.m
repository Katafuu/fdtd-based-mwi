function [epsr_final, epsr_num, epsr_den] = reconstructEpsrFinal(eps_r, angle, mask4D, maxAngle)
% reconstructEpsrFinal Pixel-wise angle-weighted TDS reconstruction.
%
% Inputs:
%   eps_r(tx,rx)        relative permittivity estimate per TX-RX pair
%   angle(tx,rx)        angle beta per TX-RX pair, radians
%   mask4D(x,y,tx,rx)   footprint/ray mask for each TX-RX pair
%   maxAngle            maximum accepted |beta|, radians
%
% Outputs:
%   epsr_final(x,y,tx)  reconstruction using each transmitter separately
%   epsr_num            numerator used for weighted averaging
%   epsr_den            denominator used for weighted averaging
%
% Correct pixel-wise formula:
%   S(x,y) = sum eps_r(tx,rx)*cos(beta)*I_txrx(x,y)
%            / sum cos(beta)*I_txrx(x,y)
%
% Note that the denominator is pixel-dependent. Using a single global
% denominator causes overlap regions to become artificially high.

    Nx_mask = size(mask4D, 1);
    Ny_mask = size(mask4D, 2);
    Ntx = size(mask4D, 3);
    Nrx = size(mask4D, 4);

    epsr_final = nan(Nx_mask, Ny_mask, Ntx);
    epsr_num   = zeros(Nx_mask, Ny_mask, Ntx);
    epsr_den   = zeros(Nx_mask, Ny_mask, Ntx);

    for tx = 1:Ntx
        numer = zeros(Nx_mask, Ny_mask);
        denom = zeros(Nx_mask, Ny_mask);

        for rx = 1:Nrx
            beta = abs(angle(tx, rx));

            if beta > maxAngle
                continue;
            end

            if ~isfinite(eps_r(tx, rx))
                continue;
            end

            I = double(mask4D(:, :, tx, rx));
            w = cos(beta);

            numer = numer + eps_r(tx, rx) * w .* I;
            denom = denom + w .* I;
        end

        epsr_tx = numer ./ denom;
        epsr_tx(denom == 0) = NaN;

        epsr_final(:, :, tx) = epsr_tx;
        epsr_num(:, :, tx) = numer;
        epsr_den(:, :, tx) = denom;
    end
end

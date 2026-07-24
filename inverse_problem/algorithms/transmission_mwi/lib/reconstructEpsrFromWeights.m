function [epsr_final, epsr_num, epsr_den] = reconstructEpsrFromWeights(eps_r, pairWeight, pairAccepted, mask4D)
% reconstructEpsrFromWeights Pixel-wise weighted TDS reconstruction.
%
% Inputs:
%   eps_r(tx,rx)          relative permittivity estimate per antenna pair
%   pairWeight(tx,rx)     scalar pair weight, e.g. cos(beta_tx)*cos(beta_rx)
%   pairAccepted(tx,rx)   logical mask selecting usable pairs
%   mask4D(x,y,tx,rx)     footprint mask for each pair
%
% Outputs:
%   epsr_final(x,y,tx)    reconstruction using each transmitter separately
%   epsr_num, epsr_den    numerator/denominator for global weighted average
%
% Formula per transmitter:
%   S_tx(x,y) = sum_rx eps_r(tx,rx)*w(tx,rx)*I_txrx(x,y)
%               / sum_rx w(tx,rx)*I_txrx(x,y)

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
            if ~pairAccepted(tx, rx)
                continue;
            end

            if ~isfinite(eps_r(tx, rx))
                continue;
            end

            w = pairWeight(tx, rx);
            if ~isfinite(w) || w <= 0
                continue;
            end

            I = double(mask4D(:,:,tx,rx));

            numer = numer + eps_r(tx, rx) * w .* I;
            denom = denom + w .* I;
        end

        epsr_tx = numer ./ denom;
        epsr_tx(denom == 0) = NaN;

        epsr_final(:,:,tx) = epsr_tx;
        epsr_num(:,:,tx) = numer;
        epsr_den(:,:,tx) = denom;
    end
end

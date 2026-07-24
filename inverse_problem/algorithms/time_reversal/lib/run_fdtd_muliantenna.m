function result = run_fdtd_muliantenna(sim, grid, material, sourceMatrix, antennas, opts)
%run_fdtd_muliantenna Run a 2-D TMz Yee FDTD simulation with all sources active.

if nargin < 6 || isempty(opts)
    opts = struct();
end

Nx = grid.sizeXY(1);
Ny = grid.sizeXY(2);
dx = grid.spacingXY(1);
dy = grid.spacingXY(2);
dt = sim.dt;
Nt = sim.Nt;
mu0 = sim.mu0;
numAntennas = antennas.numAntennas;

if ~isfield(material, 'eps_r') || ~isfield(material, 'sigma')
    error('run_fdtd_muliantenna:MissingMaterialFields', ...
        'material.eps_r and material.sigma are required.');
end

epsR = material.eps_r;
sigma = material.sigma;
if ~isfield(material, 'sigma_e') || isempty(material.sigma_e)
    sigmaE = zeros(Nx, Ny);
else
    sigmaE = material.sigma_e;
end
if ~isequal(size(epsR), [Nx Ny]) || ~isequal(size(sigma), [Nx Ny]) || ...
        ~isequal(size(sigmaE), [Nx Ny])
    error('run_fdtd_muliantenna:MaterialSizeMismatch', ...
        'material.eps_r, material.sigma, and material.sigma_e must match the grid size.');
end
if ~isnumeric(sourceMatrix) || ~isequal(size(sourceMatrix), [numAntennas Nt])
    error('run_fdtd_muliantenna:SourceSizeMismatch', ...
        'sourceMatrix must be an antennas.numAntennas-by-sim.Nt array.');
end

sigmaTotal = sigma + sigmaE;
epsMat = sim.eps0 .* epsR;
ae = (sigmaTotal .* dt) ./ (2 .* epsMat);
Ce1 = (1 - ae) ./ (1 + ae);
Ce2 = (dt ./ epsMat) ./ (1 + ae);

antennaPositions = antennas.pos;
linIdx = sub2ind([Nx Ny], antennaPositions(:, 1), antennaPositions(:, 2));
storeFieldHistory = isfield(opts, 'storeFieldHistory') && opts.storeFieldHistory;

Ez = zeros(Nx, Ny);
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
rxSignals = zeros(numAntennas, Nt);
if storeFieldHistory
    EzAll = zeros(Nx, Ny, Nt);
    HxAll = zeros(Nx, Ny, Nt);
    HyAll = zeros(Nx, Ny, Nt);
end

for n = 1:Nt
    Hx(:, 1:Ny-1) = Hx(:, 1:Ny-1) - ...
        (dt/mu0) * (Ez(:, 2:Ny) - Ez(:, 1:Ny-1)) / dy;
    Hy(1:Nx-1, :) = Hy(1:Nx-1, :) + ...
        (dt/mu0) * (Ez(2:Nx, :) - Ez(1:Nx-1, :)) / dx;

    dHyDx = (Hy(2:Nx-1, 2:Ny-1) - Hy(1:Nx-2, 2:Ny-1)) / dx;
    dHxDy = (Hx(2:Nx-1, 2:Ny-1) - Hx(2:Nx-1, 1:Ny-2)) / dy;
    curlH = dHyDx - dHxDy;
    Ez(2:Nx-1, 2:Ny-1) = ...
        Ce1(2:Nx-1, 2:Ny-1) .* Ez(2:Nx-1, 2:Ny-1) + ...
        Ce2(2:Nx-1, 2:Ny-1) .* curlH;

    for antenna = 1:numAntennas
        Ez(linIdx(antenna)) = Ez(linIdx(antenna)) + ...
            sourceMatrix(antenna, n);
    end
    rxSignals(:, n) = Ez(linIdx);

    if storeFieldHistory
        EzAll(:, :, n) = Ez;
        HxAll(:, :, n) = Hx;
        HyAll(:, :, n) = Hy;
    end
end

result = struct();
result.Ez = Ez;
result.Hx = Hx;
result.Hy = Hy;
result.rx_signals = rxSignals;
result.tx_signal = sourceMatrix;
if storeFieldHistory
    result.Ez_all = EzAll;
    result.Hx_all = HxAll;
    result.Hy_all = HyAll;
end
end

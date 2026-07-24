function plotCircularMidpoints(antenna_cell_indices, pairAccepted, Nx, Ny, pmlThickness)
% plotCircularMidpoints Plot midpoint/footprint centers of accepted pairs.

    hold on;
    Nant = size(antenna_cell_indices, 1);

    for tx = 1:Nant
        for rx = 1:Nant
            if ~pairAccepted(tx, rx)
                continue;
            end
            xip = 0.5 * (antenna_cell_indices(tx,1) + antenna_cell_indices(rx,1));
            yip = 0.5 * (antenna_cell_indices(tx,2) + antenna_cell_indices(rx,2));
            plot(xip, yip, 'k.');
        end
    end

    plot(antenna_cell_indices(:,1), antenna_cell_indices(:,2), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 3);

    rectangle('Position', [pmlThickness+1, pmlThickness+1, ...
        Nx - 2*pmlThickness - 1, Ny - 2*pmlThickness - 1], ...
        'EdgeColor', [0.7 0.7 0.7], 'LineStyle', '--');

    axis equal tight;
    xlim([1 Nx]);
    ylim([1 Ny]);
    set(gca, 'YDir', 'normal');
    xlabel('x index');
    ylabel('y index');
    grid on;
end

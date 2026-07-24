function plotAntennaArray(antenna_cell_indices, center, Nx, Ny, pmlThickness)
% plotAntennaArray Visual sanity check for antenna placement.

    plot(antenna_cell_indices(:,1), antenna_cell_indices(:,2), 'ko', ...
        'MarkerFaceColor', 'y');
    hold on;
    plot(center(1), center(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);

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

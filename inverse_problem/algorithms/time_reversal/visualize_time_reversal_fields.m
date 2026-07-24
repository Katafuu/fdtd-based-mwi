%%  Now for 2D case

Ez = E_tot;
% Ez = Ez_all;
% Ez = Ez_reversed;


nt = size(Ez, 3);

frameDelay = 0.0;
frameStride = 1;

figure;
Ez_T = permute(Ez, [2 1 3]);
EzFrame = Ez_T(:,:,1); % transpose needed for visualization only, so that x and y match on the image axis
% EzFrame = Ez(:,:,1);

h = imagesc(EzFrame);
axis equal tight;
colorbar;
colormap turbo;
clim([-0.01 0.01]);


title("E_z field, step 1");

for n = 1:frameStride:nt
    EzFrame = Ez_T(:,:,n);

    h.CData = EzFrame;

    % maxAbsFrame = max(abs(EzFrame), [], "all");
    % if maxAbsFrame == 0
    %     maxAbsFrame = 1;
    % end
    % clim([-maxAbsFrame maxAbsFrame]);

    title(sprintf("2D E_z field, step %d / %d", n, nt));
    drawnow;
    pause(frameDelay);
end


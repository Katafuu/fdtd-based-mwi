function [distance, globalAngle, beta_tx, beta_rx, pairWeight] = buildCircularPairGeometry(antenna_cell_indices, dx, dy, center)
% buildCircularPairGeometry Pairwise distances and relative antenna angles.
%
% For a circular array, each antenna is assumed to point inward toward the
% center. The angle used for weighting should therefore be an off-boresight
% angle, not the global atan2 ray angle.
%
% Outputs are Nant x Nant matrices indexed as (tx, rx):
%   distance(tx,rx)    physical distance from tx antenna to rx antenna
%   globalAngle(tx,rx) global ray angle from tx to rx
%   beta_tx(tx,rx)     outgoing ray angle relative to tx inward boresight
%   beta_rx(tx,rx)     incoming ray angle relative to rx inward boresight
%   pairWeight(tx,rx)  cos(beta_tx)*cos(beta_rx), clipped at zero

    ant_x_idx = antenna_cell_indices(:,1);
    ant_y_idx = antenna_cell_indices(:,2);

    ant_x = (ant_x_idx - 1) * dx;
    ant_y = (ant_y_idx - 1) * dy;

    center_x = (center(1) - 1) * dx;
    center_y = (center(2) - 1) * dy;

    delta_x = ant_x.' - ant_x;  % tx row, rx column
    delta_y = ant_y.' - ant_y;

    distance = sqrt(delta_x.^2 + delta_y.^2);
    globalAngle = atan2(delta_y, delta_x);

    % Inward boresight angle of each antenna: antenna -> center.
    boresight = atan2(center_y - ant_y, center_x - ant_x);  % Nant x 1

    % TX side: outgoing direction tx -> rx compared to tx boresight.
    beta_tx = localWrapToPi(globalAngle - boresight);

    % RX side: incoming ray is viewed as rx -> tx, so reverse the direction.
    reverseAngle = atan2(-delta_y, -delta_x);
    beta_rx = localWrapToPi(reverseAngle - boresight.');

    pairWeight = max(cos(beta_tx), 0) .* max(cos(beta_rx), 0);

    d = eye(size(distance)) == 1;
    distance(d) = NaN;
    globalAngle(d) = NaN;
    beta_tx(d) = NaN;
    beta_rx(d) = NaN;
    pairWeight(d) = 0;
end

function a = localWrapToPi(a)
% localWrapToPi Toolbox-free wrap to [-pi, pi].
    a = atan2(sin(a), cos(a));
end

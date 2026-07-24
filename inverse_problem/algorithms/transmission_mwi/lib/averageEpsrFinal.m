function epsr_avg = averageEpsrFinal(epsr_num, epsr_den)
% averageEpsrFinal Correct global weighted average over transmitters.
%
% This is preferable to mean(epsr_final,3), because it preserves the
% pixel-wise denominators/weights from every transmitter-receiver pair.

    den_total = sum(epsr_den, 3);
    num_total = sum(epsr_num, 3);

    epsr_avg = num_total ./ den_total;
    epsr_avg(den_total == 0) = NaN;
end

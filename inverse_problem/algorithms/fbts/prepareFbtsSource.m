function [cfg, time, pulse, weight, source] = prepareFbtsSource(cfg)
%prepareFbtsSource Single definition of the FBTS excitation and time weight.
validateattributes(cfg.Nt, {'numeric'}, {'scalar','integer','>=',2});
validateattributes(cfg.dt, {'numeric'}, {'scalar','positive','finite'});
time = (0:cfg.Nt-1) .* cfg.dt;
tau = 0.125e-9;
cfg.source.func = @(t) ...
    (4 .* t.^3 ./ tau.^4 - t.^4 ./ tau.^5) .* exp(-t ./ tau);
pulse = cfg.source.func(time);
weight = cos(pi .* time ./ (2 .* time(end)));
weight(end) = 0;
source = struct('name', 'paper_pulse', 'tau_seconds', tau, ...
    'formula', '(4*t^3/tau^4-t^4/tau^5)*exp(-t/tau)', ...
    'injection', 'additive Ez point source', 'signal', 'total Ez', ...
    'amplitude_units', 'solver Ez units; no experimental calibration', ...
    'self_receiver_included', true);
end

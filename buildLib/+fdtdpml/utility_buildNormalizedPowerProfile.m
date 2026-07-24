function profile = utility_buildNormalizedPowerProfile(depth, thickness, order)
%utility_buildNormalizedPowerProfile Compute a clamped power-law PML profile.

if ~isnumeric(depth) || any(~isfinite(depth(:)))
    error('fdtdpml:utility_buildNormalizedPowerProfile:InvalidDepth', 'depth must be a finite numeric array.');
end
if ~isscalar(thickness) || thickness <= 0 || ~isfinite(thickness)
    error('fdtdpml:utility_buildNormalizedPowerProfile:InvalidThickness', 'thickness must be a positive finite scalar.');
end
if ~isscalar(order) || order <= 0 || ~isfinite(order)
    error('fdtdpml:utility_buildNormalizedPowerProfile:InvalidOrder', 'order must be a positive finite scalar.');
end

normalizedDepth = min(max(depth, 0), thickness) / thickness;
profile = normalizedDepth .^ order;
end

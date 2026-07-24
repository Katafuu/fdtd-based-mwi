function sigmaMax = utility_sigmaMaxFromReflection(order, targetReflection, eta0, physicalThickness)
%utility_sigmaMaxFromReflection Compute the maximum PML conductivity.

if ~isscalar(order) || order <= 0 || ~isfinite(order)
    error('fdtdpml:utility_sigmaMaxFromReflection:InvalidOrder', 'order must be a positive finite scalar.');
end
if ~isscalar(targetReflection) || targetReflection <= 0 || targetReflection >= 1
    error('fdtdpml:utility_sigmaMaxFromReflection:InvalidReflection', ...
        'targetReflection must be between 0 and 1.');
end
if ~isscalar(eta0) || eta0 <= 0 || ~isfinite(eta0)
    error('fdtdpml:utility_sigmaMaxFromReflection:InvalidEta0', 'eta0 must be a positive finite scalar.');
end
if ~isscalar(physicalThickness) || physicalThickness <= 0 || ~isfinite(physicalThickness)
    error('fdtdpml:utility_sigmaMaxFromReflection:InvalidThickness', ...
        'physicalThickness must be a positive finite scalar.');
end

sigmaMax = -(order + 1) * log(targetReflection) / (2 * eta0 * physicalThickness);
end

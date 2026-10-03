function reconstruction = reconstructHomogeneousTargetMap( ...
        targetMask, backgroundEpsr, recoveredTargetEpsr)
%reconstructHomogeneousTargetMap Fill a known target mask with one value.

validateattributes(targetMask, {'logical'}, {'2d', 'nonempty'}, ...
    mfilename, 'targetMask');
validateattributes(backgroundEpsr, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'backgroundEpsr');
validateattributes(recoveredTargetEpsr, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'recoveredTargetEpsr');

reconstruction = backgroundEpsr .* ones(size(targetMask));
reconstruction(targetMask) = recoveredTargetEpsr;
end

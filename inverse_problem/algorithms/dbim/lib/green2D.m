function greenMatrix = green2D(observationPositions, sourcePositions, wavenumber)
%green2D Compute the outgoing 2-D scalar Green function.
%   greenMatrix(i,j) is the Green function from sourcePositions(j,:) to
%   observationPositions(i,:), using the exp(+1i*omega*t) convention.

observationPositions = validatePositions(observationPositions, ...
    'green2D:InvalidObservationPositions');
sourcePositions = validatePositions(sourcePositions, ...
    'green2D:InvalidSourcePositions');
if ~isnumeric(wavenumber) || ~isscalar(wavenumber) || ...
        ~isfinite(wavenumber) || wavenumber == 0
    error('green2D:InvalidWavenumber', ...
        'wavenumber must be a finite nonzero numeric scalar.');
end
wavenumber = double(wavenumber);

deltaX = observationPositions(:, 1) - sourcePositions(:, 1).';
deltaY = observationPositions(:, 2) - sourcePositions(:, 2).';
distance = hypot(deltaX, deltaY);

% Preserve a finite cell-centre approximation if source and observation
% coordinates coincide. DBIM normally keeps antennas outside the DOI.
coincident = distance == 0;
distance(coincident) = eps / max(1, abs(wavenumber));

greenMatrix = (-1i/4) * besselh(0, 2, wavenumber * distance);
if any(~isfinite(greenMatrix(:)))
    error('green2D:NonfiniteResult', ...
        'The Green function evaluation produced a nonfinite value.');
end
end

function positions = validatePositions(positions, errorId)
if ~isnumeric(positions) || size(positions, 2) ~= 2 || ...
        isempty(positions) || any(~isfinite(positions(:)))
    error(errorId, 'Positions must be a nonempty finite N-by-2 array.');
end
positions = double(positions);
end

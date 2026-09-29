% gen_combinatorial_scenarios Build target positions and antenna-off layouts.
% Run build_cfg.m first. Add the desired physical [x y] positions in meters below.
targetPositions = [
    0.146000  0.138000;
    0.152000  0.138000;
    0.158000  0.138000;
    0.164000  0.138000;
    0.170000  0.138000;
    0.176000  0.138000;
    0.182000  0.138000;
    0.188000  0.138000;
    0.194000  0.138000;
    0.200000  0.138000;
    0.200000  0.132000;
    0.194000  0.132000;
    0.188000  0.132000;
    0.182000  0.132000;
];
antennaOffConfigs = generateAntennaOffConfigs(cfg.antennas.numAntennas);

function configs = generateAntennaOffConfigs(numAntennas)
% Generate binary bracelets without materializing all 2^N antenna subsets.
configs = cell(1, numAntennas);
configs{1} = zeros(1, 0); % One all-on configuration.
for offCount = 1:numAntennas-1
    configs{offCount + 1} = zeros(0, offCount);
end

necklace = zeros(1, numAntennas + 1);
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
enumerateNecklaces(1, 1);
for offCount = 1:numAntennas-1
    configs{offCount + 1} = sortrows(configs{offCount + 1});
end

    function enumerateNecklaces(position, period)
        if position > numAntennas
            if mod(numAntennas, period) == 0
                retainBracelet(necklace(2:end));
            end
            return
        end

        necklace(position + 1) = necklace(position - period + 1);
        enumerateNecklaces(position + 1, period);
        if necklace(position - period + 1) == 0
            necklace(position + 1) = 1;
            enumerateNecklaces(position + 1, position);
        end
    end

    function retainBracelet(disabledBits)
        offCount = sum(disabledBits);
        if offCount == 0 || offCount == numAntennas
            return
        end

        canonicalKey = '';
        representative = [];
        for shift = 0:numAntennas-1
            for reflected = [false true]
                if reflected
                    order = mod(shift - (0:numAntennas-1), numAntennas) + 1;
                else
                    order = mod(shift + (0:numAntennas-1), numAntennas) + 1;
                end
                transformed = disabledBits(order);
                candidateKey = char(double('0') + transformed);
                if isempty(canonicalKey)
                    canonicalKey = candidateKey;
                else
                    orderedKeys = sort({canonicalKey, candidateKey});
                    canonicalKey = orderedKeys{1};
                end

                candidateIndices = find(transformed);
                if isempty(representative)
                    representative = candidateIndices;
                else
                    firstDifference = find(candidateIndices ~= representative, 1);
                    if ~isempty(firstDifference) && ...
                            candidateIndices(firstDifference) < ...
                            representative(firstDifference)
                        representative = candidateIndices;
                    end
                end
            end
        end

        if ~isKey(seen, canonicalKey)
            seen(canonicalKey) = true;
            configs{offCount + 1}(end + 1, :) = representative;
        end
    end
end

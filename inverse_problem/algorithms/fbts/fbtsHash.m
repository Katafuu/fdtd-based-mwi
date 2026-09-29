function digest = fbtsHash(value, mode)
%fbtsHash SHA-256 of typed values (sorted struct fields), or raw file bytes.
md = java.security.MessageDigest.getInstance('SHA-256');
if nargin > 1 && strcmp(mode, 'file')
    fid = fopen(value, 'rb');
    assert(fid >= 0, 'fbts:HashFile', 'Cannot read %s.', value);
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    while ~feof(fid)
        chunk = fread(fid, 1048576, '*uint8');
        if ~isempty(chunk), md.update(typecast(chunk(:), 'int8')); end
    end
else
    append(value);
end
digest = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));

    function bytes(data)
        if ~isempty(data), md.update(typecast(uint8(data(:)), 'int8')); end
    end
    function append(v)
        bytes(unicode2native([class(v) ':'], 'UTF-8'));
        bytes(typecast(uint64([numel(size(v)), size(v)]), 'uint8'));
        if isstruct(v)
            names = sort(fieldnames(v));
            append(names);
            for element = 1:numel(v)
                for field = 1:numel(names)
                    append(v(element).(names{field}));
                end
            end
        elseif iscell(v)
            for element = 1:numel(v), append(v{element}); end
        elseif isstring(v)
            assert(~any(ismissing(v(:))), 'fbts:HashType', 'Missing strings cannot be hashed.');
            append(cellstr(v));
        elseif ischar(v)
            bytes(typecast(uint16(v(:)), 'uint8'));
        elseif isnumeric(v)
            assert(~issparse(v), 'fbts:HashType', 'Sparse values are not supported.');
            bytes(uint8(~isreal(v)));
            bytes(typecast(real(v(:)), 'uint8'));
            if ~isreal(v), bytes(typecast(imag(v(:)), 'uint8')); end
        elseif islogical(v)
            bytes(uint8(v));
        else
            error('fbts:HashType', 'Unsupported fingerprint type: %s.', class(v));
        end
    end
end

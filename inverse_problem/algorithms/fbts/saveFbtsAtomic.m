function saveFbtsAtomic(path, variables)
%saveFbtsAtomic Commit top-level MAT v7.3 variables using a same-folder rename.
folder = fileparts(path);
if ~isfolder(folder), mkdir(folder); end
temporary = [tempname(folder) '.mat'];
cleanup = onCleanup(@() removeTemporary(temporary)); %#ok<NASGU>
save(temporary, '-struct', 'variables', '-v7.3');
[ok,message] = movefile(temporary,path,'f');
assert(ok,'fbts:CommitFailed','%s',message);
end
function removeTemporary(path)
if isfile(path), delete(path); end
end

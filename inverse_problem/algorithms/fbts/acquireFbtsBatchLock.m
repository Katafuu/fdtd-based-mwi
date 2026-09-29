function cleanup = acquireFbtsBatchLock(root)
%acquireFbtsBatchLock OS-backed lock released automatically on process death.
file = java.io.RandomAccessFile(fullfile(root,'.batch.lock'),'rw');
channel = file.getChannel();
try
    lock = channel.tryLock();
catch
    channel.close(); file.close();
    error('fbts:BatchLocked','Another coordinator is using this batch directory.');
end
if isempty(lock)
    channel.close(); file.close();
    error('fbts:BatchLocked','Another coordinator is using this batch directory.');
end
cleanup = onCleanup(@release);
    function release()
        lock.release(); channel.close(); file.close();
    end
end

function runFbtsQueue(tasks,shared,width,onComplete)
%runFbtsQueue Rolling bounded queue; the coordinator is the only index writer.
if width == 1
    for k = 1:numel(tasks)
        tasks{k}.submitted = datetime('now','TimeZone','UTC');
        onComplete(k,runFbtsQueuedCase(tasks{k},shared));
    end
    return
end
constant = parallel.pool.Constant(shared);
futures = parallel.FevalFuture.empty; indices = []; next = 1;
cleanup = onCleanup(@release); %#ok<NASGU>
while next <= numel(tasks) || ~isempty(futures)
    while next <= numel(tasks) && numel(futures) < width
        tasks{next}.submitted = datetime('now','TimeZone','UTC');
        futures(end+1) = parfeval(@runFbtsQueuedCase,1,tasks{next},constant);
        indices(end+1) = next; next = next+1;
    end
    try
        [finished,row] = fetchNext(futures);
    catch exception
        finished = find(arrayfun(@(f) strcmp(f.State,'finished') && ~isempty(f.Error),futures),1);
        if isempty(finished), rethrow(exception); end
        row = tasks{indices(finished)}.row;
        row.status = "failed"; row.failure_stage = "worker";
        row.error_id = string(exception.identifier); row.error_message = string(exception.message);
    end
    onComplete(indices(finished),row);
    futures(finished) = []; indices(finished) = [];
end
    function release()
        if ~isempty(futures), cancel(futures); end
        delete(constant);
    end
end

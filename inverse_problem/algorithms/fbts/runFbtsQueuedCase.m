function row = runFbtsQueuedCase(task,shared)
%runFbtsQueuedCase Convert ordinary failures into results for fetchNext.
if isa(shared,'parallel.pool.Constant'), shared = shared.Value; end
task.cfg = shared.cfg; task.measurements = shared.measurements;
task.row.queue_wait_seconds = seconds(datetime('now','TimeZone','UTC')-task.submitted);
try
    row = runFbtsBatchCase(task);
catch exception
    row = task.row; row.status = "failed"; row.failure_stage = "worker_or_status_write";
    row.error_id = string(exception.identifier); row.error_message = string(exception.message);
end
end

function value = fbtsUtcNow()
%fbtsUtcNow UTC chronology, separate from elapsed-duration timers.
value=string(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

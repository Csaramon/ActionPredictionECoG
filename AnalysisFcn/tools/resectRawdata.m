time2remove = 83;  % unit in second
time2removeEnd = 90;  % unit in second
timept = raw_data.fsample*time2remove;
timeptEnd = raw_data.fsample*time2removeEnd;

raw_data.trial{1,1}(:,1:timept) = [];
raw_data.time{1,1}(:,end-timept+1:end) = [];
raw_data.trial{1,1}(:,end-timeptEnd+1:end) = [];
raw_data.time{1,1}(:,end-timeptEnd+1:end) = [];

raw_data.sampleinfo(2) = size(raw_data.trial{1,1},2);
raw_data.cfg.trl(2) = size(raw_data.trial{1,1},2);
raw_data.hdr.nSamples = size(raw_data.trial{1,1},2);

eventdata.latency = eventdata.latency-time2remove*1000;
eventdata.orig_latency = eventdata.orig_latency-time2remove*1000;

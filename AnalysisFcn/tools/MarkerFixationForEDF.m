%% setction 0
resortInd = [];
for il = 1:numel(lijunwen_label)
    tmpInd = find(strcmp(raw_data.label,lijunwen_label{il}));
    
    resortInd = [resortInd,tmpInd];
    
    
end

raw_data.label = lijunwen_label;
 raw_data.hdr.nChans = numel(lijunwen_label);
        raw_data.hdr.label = lijunwen_label;
        raw_data.hdr.chanunit = raw_data.hdr.chanunit(resortInd);
        raw_data.hdr.chantype = raw_data.hdr.chantype(resortInd);
        raw_data.cfg.channel = lijunwen_label;
raw_data.trial{1} = raw_data.trial{1}(resortInd,:);



%% setction 1
minMarkInter = 0.12; % minimum button press interval in seconds

[befile, bepath, bestatus] = uigetfile({'*.mat';},'Please choose the behavioral data.');

if  bestatus==1
    
    aa = load([bepath befile]);
    Behavior = aa.RECORD;
    
    true_marker = Behavior.answer(1);
    true_eventinfo = Behavior.timeStamp(1);
    if ~isfield(Behavior,'ITrial')
        Behavior.ITrial = Behavior.iTrial;
    end
    true_trial = Behavior.ITrial(1);
    
    
    
    for itrl = unique(Behavior.ITrial)'
        trlInd = find(Behavior.ITrial==itrl);
        
        trialAns = Behavior.answer(trlInd);
        trialTime = Behavior.timeStamp(trlInd);
        trialInd = Behavior.ITrial(trlInd);
        
        n =2;
        while n <= numel(trialAns)
            if (trialAns(n)==trialAns(n-1)) & ...
                    (trialTime(n)-trialTime(n-1)) > minMarkInter
                true_marker = [true_marker;trialAns(n)];
                true_eventinfo = [true_eventinfo;trialTime(n)];
                true_trial = [true_trial;trialInd(n)];
            end
            if (trialAns(n)~=trialAns(n-1))
                true_marker = [true_marker;trialAns(n)];
                true_eventinfo = [true_eventinfo;trialTime(n)];
                true_trial = [true_trial;trialInd(n)];
            end
            n=n+1;
        end
        
        
        
    end
    
end
toinspectbe = [true_eventinfo,true_marker,true_trial];

%% setction 2

toinspectEEG(:,1) = mat2cell(eventdata.marker,ones(size(eventdata.marker,1),1),size(eventdata.marker,2));
toinspectEEG(:,2) = mat2cell(eventdata.latency./1000,ones(size(eventdata.latency,1),1),size(eventdata.latency,2));

% replay condition
toinspectbe(toinspectbe(:,2)==7,:) = [];

% % part 1
toinspectbe1 = toinspectbe(toinspectbe(:,3)==1,:);
toinspectbe2 = toinspectbe(toinspectbe(:,3)==2,:);
toinspectbe3 = toinspectbe(toinspectbe(:,3)==3,:);
toinspectbe4 = toinspectbe(toinspectbe(:,3)==4,:);
toinspectbe5 = toinspectbe(toinspectbe(:,3)==5,:);
toinspectbe6 = toinspectbe(toinspectbe(:,3)==6,:);



% % part 2
% need user action
toinspectbe1(:,1) = toinspectbe1(:,1)+toinspectEEG{3,2}-toinspectbe1(1,1);
% need user action
toinspectbe2(:,1) = toinspectbe2(:,1)+toinspectEEG{16,2}-toinspectbe2(1,1);
% need user action
toinspectbe3(:,1) = toinspectbe3(:,1)+toinspectEEG{34,2}-toinspectbe3(1,1);
% need user action
toinspectbe4(:,1) = toinspectbe4(:,1)+toinspectEEG{48,2}-toinspectbe4(1,1);
% need user action
toinspectbe5(:,1) = toinspectbe5(:,1)+toinspectEEG{83,2}-toinspectbe5(1,1);
% need user action
toinspectbe6(:,1) = toinspectbe6(:,1)+toinspectEEG{122,2}-toinspectbe6(1,1);

orig_marker = cell2mat(toinspectEEG(:,1));
orig_latency = cell2mat(toinspectEEG(:,2)).*1000;



if size(orig_latency,1) < size(orig_latency,2)
    orig_latency = orig_latency';
end

eventdata.orig_marker = orig_marker;
eventdata.orig_latency = orig_latency;
eventdata.marker = orig_marker;
eventdata.latency = orig_latency;

uisave('eventdata', [bepath '_eventdata.mat']);








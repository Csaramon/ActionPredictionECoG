cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [2 6];
cfg.bpfilttype = 'firws';

filtRerefData = ft_preprocessing(cfg,rerefData);

eeglab;close
%%%%%%%%%%%%%%%%%%%
channel2check = [38,39,40,41,43,47];
data2check = filtRerefData.trial{1}(channel2check,:);



for i = 1:numel(eventdata.latency)
    
    EEG.event(i).type = eventdata.marker(i,:);
    
    EEG.event(i).latency = eventdata.latency(i)*2;
      
end


eegplot(data2check,'srate',2000,'events',EEG.event);

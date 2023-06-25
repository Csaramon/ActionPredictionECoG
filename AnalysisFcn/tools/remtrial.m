function badtrialindex = remtrial(epochdata,thr) 

low = nanmean(epochdata,3)-thr*nanstd(epochdata,0,3); % the zscore across trials; 2D data: electrodes*timepoints
high = nanmean(epochdata,3)+thr*nanstd(epochdata,0,3);
outlier1 = (epochdata < repmat(low,[1 1 size(epochdata,3)]));
outlier2 = (epochdata > repmat(high,[1 1 size(epochdata,3)]));
outlier = outlier1 | outlier2;  % electrodes*timepoint*trials
badtrialindex = find(sum(squeeze(sum(outlier,1)),1)>(0.01*size(epochdata,2)));
% first squeeze the electrodes and then squeeze timepoint,finally get the badindex in trial dimension
% the sum(X,D) function will decrease the dimension D to 1, along with the squeeze function will remove the dimension D, that is get the sum matrix(without dimension D) along the dimension D

% to display the badtrial index in each electrode
% chan_n_trial = squeeze(sum(outlier,2)); %=% electrode*trial; number of badtrials in each electrode
%  for i =  1:size(epochdata,1) %=% electrode
%      if sum(chan_n_trial(i,:))>0
%      badtrial_index = find(chan_n_trial(i,:));
%      disp(['badtrials for elec_' num2str(i) ':' num2str(badtrial_index)]);
%      end
%  end
% end
%  
%  
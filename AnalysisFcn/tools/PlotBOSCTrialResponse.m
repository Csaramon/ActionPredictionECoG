function PlotBOSCTrialResponse


latency = [-1000,0]; % 2s at 500Hz sampling rate
foi = [2 8];

 [FILENAME, PATHNAME, FILTERINDEX] = uigetfile(pwd, 'Pick a BOSC data file');
 
 if FILTERINDEX >= 0
            tempData=  load([PATHNAME FILENAME]);
            a = fieldnames(tempData);
            BOSC = tempData.(a{1});
 end

 Findex = BOSC.F_b>=foi(1) & BOSC.F_b<=foi(2);
 normaliseData = zscore(BOSC.B_b(Findex,:),0,2);
 threshData = BOSC.B_b(Findex,:).*BOSC.detected_b(Findex,:);
%  normaliseData = mean((BOSC.B_b(Findex,:)-repmat(BOSC.Bmean_b(Findex),1,size(BOSC.B_b,2))) ...
%      ./repmat(BOSC.Bmean_b(Findex),1,size(BOSC.B_b,2)),1);
 trlLatency = round(BOSC.eventdata_b.orig_latency./2);
 trlLatency(trlLatency < abs(latency(1)) | trlLatency > abs(latency(2)) + size(BOSC.B_b,2)) = [];
 
  trlDataAll = [];
 for itrl = 1:size(trlLatency)
     
     trlData = normaliseData(:,trlLatency(itrl) + latency(1):trlLatency(itrl) + latency(2));
     thrCrit = threshData(:,trlLatency(itrl) + latency(1):trlLatency(itrl) + latency(2));
     trlData(mean(trlData,2)==0,:) = [];
     trlDataFilt = mean(trlData,1);
     if max(trlDataFilt) > 1.96
     trlDataAll = [trlDataAll;trlDataFilt];
     
     end
 end
  
 [maxResp,maxInd] = max(trlDataAll,[],2);
[~,IA] = sort(maxInd);

trlDataAll = trlDataAll(IA,:);

 subplot(4,1,1:3)
 imagesc(trlDataAll)
 caxis([-1 3])
 ylabel('Trial Number')
 set(gca,'xtick',[])
  title('Trial response of sample electrode')
 subplot(4,1,4)
 
 plot(-2:1/500:0,sum(trlDataAll,1)./size(trlDataAll,1),'k','linewidth',1.5)
 set(gca, 'TickDir','out')
 set(gca, 'LineWidth',1.5)
 xlabel('Time relative to response (s)')
 box off
 
 
 
end

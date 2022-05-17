allMetric = cat(1,allMetricM,allMetricS);
allSubInd = unique(lmeTBL.Sub);

toi = Para.timePT > 0.45 & Para.timePT < 0.55;
x2plot = Para.freq;
for isub = 1:numel(allSubInd)
    
    elecIndM  = find(lmeTBL.Sub==allSubInd(isub) & lmeTBL.Cond=='1');
    elecIndS  = find(lmeTBL.Sub==allSubInd(isub) & lmeTBL.Cond=='2');
    
MetricM = squeeze(nanmean(allMetric(elecIndM,:,toi),3));
MetricS = squeeze(nanmean(allMetric(elecIndS,:,toi),3));

figure;

% subplot(1,2,1)
% hold on
% plot(MetricM')
% plot(MetricS')
% subplot(1,2,2)
hold on
hM = shadedErrorBar(x2plot, mean(MetricM),std(MetricM)./sqrt(size(MetricM,1)),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(x2plot, mean(MetricS), std(MetricS)./sqrt(size(MetricS,1)),{'color',[30 144 255]/255},1);
ylabel('Coherence (a.u.)')
xlabel('Frequency (Hz)')
xlim([0 120])
    
% saveas(gcf,[pathname 'sub' num2str(isub) 'COH'])
end
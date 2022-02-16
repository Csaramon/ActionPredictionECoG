%% plot time-frequency data after reducing dimension


% foi  = Para.freq>60 & Para.freq<90;
% foi  = Para.freq>20 & Para.freq<30;
% MetricM = squeeze(mean(allMetricM(:,foi,:),2));
% MetricS = squeeze(mean(allMetricS(:,foi,:),2));
% x2plot = Para.timePT;

toi = Para.timePT > 0 & Para.timePT < 0.5;
MetricM = squeeze(nanmean(allMetricM(:,:,toi),3));
MetricS = squeeze(nanmean(allMetricS(:,:,toi),3));
x2plot = Para.freq;


tMap = [];
pMap = [];
y2plot = [];
se2plot = [];
yraw = [];
seraw = [];


strlen = 0;

for itime = 1:size(MetricM,2)
    
    s = ['Calculating tf point:' num2str(itime) '/' num2str(size(MetricM,2)) 'points'];
    strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
    strlen = strlentmp - strlen;
    
    frameDataM = double(MetricM(:,itime));
    frameDataS = double(MetricS(:,itime));
    % skip nan point
    if ~any(frameDataM,'all')
        continue
    end
    
    lmeTBL.Y = [frameDataM;frameDataS];
    
    
    lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
    
%     lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec1)+(1|Elec2)+(1|Elec1:Elec2)','fitmethod','reml','DummyVarCoding','effects');
    
    [~,~,lmeStats] = fixedEffects(lmeStruct);
    tMap(1,itime) = lmeStats.tStat(2);
    pMap(1,itime) = lmeStats.pValue(2);
    [randBeta,~,~] = randomEffects(lmeStruct);
    Z = designMatrix(lmeStruct,'random');
    Ycorr = lmeTBL.Y-Z*randBeta;
    obsVal1 = Ycorr(lmeTBL.Cond=='1');
    obsVal2 = Ycorr(lmeTBL.Cond=='2');
    y2plot(:,itime) = [mean(obsVal1);mean(obsVal2)];
    se2plot(:,itime) = [std(obsVal1)./sqrt(numel(obsVal1));std(obsVal2)./sqrt(numel(obsVal2))];
    rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
    rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
    yraw(:,itime) = [mean(rawVal1);mean(rawVal2)];
    seraw(:,itime) = [std(rawVal1)./sqrt(numel(rawVal1));std(rawVal2)./sqrt(numel(rawVal2))];
    
end



% plot new result

y2plot = yraw;
se2plot = seraw;

x2plot = x2plot;

% FDR  correction
if numel(x2plot) == numel(Para.timePT)
    corrTimeInd = Para.timePT>-0.55 & Para.timePT< 1.05;
    [p_fdr, p_masked] = fdr(pMap(corrTimeInd), 0.05,'Parametric');
    highlight = pMap<=p_fdr;
else
    corrFreqInd = Para.freq>2 & Para.freq< 120;
    [p_fdr, p_masked] = fdr(pMap(corrFreqInd), 0.05,'Parametric');
    highlight = pMap<=p_fdr;
end

highlight =double(highlight);
highlight(highlight==0) = nan;

hf = figure;

hold on
hM = shadedErrorBar(x2plot, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(x2plot, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);

hsig = plot(x2plot,(max(y2plot(:))+0.2*range(y2plot(:)))*highlight,'k-','linewidth',2);
title(filename(1:end-4))
% xlim([-0.5 1])
% ylim([-2*max(abs(y2plot(:))),2*max(abs(y2plot(:)))])
plot([0,0],get(gca,'ylim'),'k--')
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);

ti = find(filename=='_');
xlabel('Time relative to camera change (sec)')
ylabel([filename(1:ti-1) ' from    PSI    from ' filename(ti+1:end-4)])



      

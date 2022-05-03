%% plot ERP

clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/BP/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap);

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;


% the significant points could be marked with a red stars
highlight =double(highlight);
highlight(highlight==0) = nan;

hf = figure;

hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);

hsig = plot(Para.time,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'k*');

title(filename(1:end-4))
xlabel('Time relative to camera change (sec)')
ylabel('ERP (μV)')
xlim([-0.5 1])
ylim([ -3 3])
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')

legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);

% save the figure to data location
saveas(hf,[pathname filename(1:end-4)])

%% plot Power Spectrum
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/BP/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(tMap==0)=nan;

lowFind = Para.freq <= 30;
highFind = Para.freq > 30;
% normalise
y2plot = y2plot.*Para.freq;
se2plot = se2plot.*Para.freq;

% uncorrected
highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap);

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlightcorr = p_masked;

% the significant points could be marked with a red stars
highlight =double(highlight);
highlightcorr =double(highlightcorr);
highlight(highlight==0) = nan;
highlightcorr(highlightcorr==0) = nan;

hf = figure('Position', [400 400 800 600]);
set(gca,'linewidth',2)
hold on
xlim([0 120])


yyaxis left
set(gca,'ycolor','k')
ylabel('Low Frequency Power (a.u.)')
hlowM = shadedErrorBar(Para.freq(lowFind), y2plot(1,lowFind),se2plot(1,lowFind),{'color',[255 106 106]/255,'linewidth',1,'linestyle','-'},1);
hlowS = shadedErrorBar(Para.freq(lowFind), y2plot(2,lowFind), se2plot(2,lowFind),{'color',[30 144 255]/255,'linewidth',1,'linestyle','-'},1);
hlowS.mainLine.Marker = 'none';

hlowsig = plot(Para.freq(lowFind),(max(y2plot(:,lowFind),[],'all')+0.2*range(y2plot(:,lowFind),'all'))*highlight(lowFind),'-','color',[0.5 0.5 0.5],'LineWidth',2);
hlowsigcorr = plot(Para.freq(lowFind),(max(y2plot(:,lowFind),[],'all')+0.1*range(y2plot(:,lowFind),'all'))*highlightcorr(lowFind),'-','color',[0 0 0],'LineWidth',2);

yyaxis right
set(gca,'ycolor','k')
ylabel('High Frequency Power (a.u.)')
hhighM = shadedErrorBar(Para.freq(highFind), y2plot(1,highFind),se2plot(1,highFind),{'color',[255 106 106]/255,'linewidth',1,'linestyle','-'},1);
hhighS = shadedErrorBar(Para.freq(highFind), y2plot(2,highFind), se2plot(2,highFind),{'color',[30 144 255]/255,'linewidth',1,'linestyle','-'},1);
hhighS.mainLine.Marker = 'none';

hhighsig = plot(Para.freq(highFind),(max(y2plot(:,highFind),[],'all')+0.2*range(y2plot(:,highFind),'all'))*highlight(highFind),'-','color',[0.5 0.5 0.5],'LineWidth',2);
hhighsigcorr = plot(Para.freq(highFind),(max(y2plot(:,highFind),[],'all')+0.1*range(y2plot(:,highFind),'all'))*highlightcorr(highFind),'-','color',[0 0 0],'LineWidth',2);


title(filename(1:end-4))
xlabel('Frequency (Hz)')

legend([hlowM.mainLine,hlowS.mainLine,hlowsig,hlowsigcorr], ...
    ['Intact'],['Scrambled'],['P<0.05'],['P<0.05 (corrected)'],'box','off','NumColumns',2);


% plot the inlets of beta and gamma power of each participant
subColor = distinguishable_colors(20);
subColor = subColor([1 2 3 5 7 9 11 15 16 18],:);

betaF = [20 30];
gammaF = [60 90];
% betaInd = Para.freq>=min(betaF)-0.5 &  Para.freq<=max(betaF)+0.5;
% gammaInd = Para.freq>=min(gammaF)-0.5 &  Para.freq<=max(gammaF)+0.5;
betaInd = Para.freq>min(betaF) &  Para.freq<max(betaF);
gammaInd = Para.freq>min(gammaF) &  Para.freq<max(gammaF);

betaPowM = mean(allMetricM(:,betaInd),2);
betaPowS = mean(allMetricS(:,betaInd),2);

gammaPowM = mean(allMetricM(:,gammaInd),2);
gammaPowS = mean(allMetricS(:,gammaInd),2);

haBeta = axes(hf,'Position',[0.5 0.5 0.15 0.25],...
    'Xlim',[0.5 2.5],'XTick', [1,2],'XTickLabel',{'I','S'},'YTick',[],'NextPlot','add');
ylabel(haBeta,'Beta Power')
haGamma = axes(hf,'Position',[0.7 0.5 0.15 0.25],...
    'Xlim',[0.5 2.5],'XTick', [1,2],'XTickLabel',{'I','S'},'YTick',[],'NextPlot','add');
ylabel(haGamma,'Gamma Power')


% compare the mean value of each frequecy band
lmeTBL.Y(lmeTBL.Cond=='1') = betaPowM;
lmeTBL.Y(lmeTBL.Cond=='2') = betaPowS;
lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
[~,~,lmeStatsBeta] = fixedEffects(lmeStruct);
% subtract the random effect
[randBeta,~,~] = randomEffects(lmeStruct);
Z = designMatrix(lmeStruct,'random');
YcorrBeta = lmeTBL.Y-Z*randBeta;

lmeTBL.Y(lmeTBL.Cond=='1') =gammaPowM;
lmeTBL.Y(lmeTBL.Cond=='2') = gammaPowS;
lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
[~,~,lmeStatsGamma] = fixedEffects(lmeStruct);
% subtract the random effect
[randBeta,~,~] = randomEffects(lmeStruct);
Z = designMatrix(lmeStruct,'random');
YcorrGamma = lmeTBL.Y-Z*randBeta;

ysubbetaM = [];
ysubbetaS = [];
ysubgammaM = [];
ysubgammaS = [];
for isub = unique(lmeTBL.Sub)'
    
    ybetaM = mean(YcorrBeta(lmeTBL.Cond=='1' & lmeTBL.Sub==isub));
    ybetaS = mean(YcorrBeta(lmeTBL.Cond=='2' & lmeTBL.Sub==isub));
    plot(haBeta,[1,2],[ybetaM,ybetaS],'-','Color',[subColor(double(isub),:) 0.3])
    
    ygammaM = mean(YcorrGamma(lmeTBL.Cond=='1' & lmeTBL.Sub==isub));
    ygammaS = mean(YcorrGamma(lmeTBL.Cond=='2' & lmeTBL.Sub==isub));
    plot(haGamma,[1,2],[ygammaM,ygammaS],'-','Color',[subColor(double(isub),:) 0.3])
    
    ysubbetaM = [ysubbetaM,ybetaM];
    ysubbetaS = [ysubbetaS,ybetaS];
    ysubgammaM = [ysubgammaM,ygammaM];
    ysubgammaS = [ysubgammaS,ygammaS];
end

% plot mean value
plot(haBeta,[1,2],[nanmean(ysubbetaM),nanmean(ysubbetaS)],'-','Color',[0.4 0.4 0.4 1],'LineWidth',2)
set(haBeta,'LineWidth',1.5)
% a = get(haBeta,'Children');
% set(haBeta,'Children',[a(2:end);a(1)])
plot(haGamma,[1,2],[nanmean(ysubgammaM),nanmean(ysubgammaS)],'-','Color',[0.4 0.4 0.4 1],'LineWidth',2)
set(haGamma,'LineWidth',1.5)
% a = get(haGamma,'Children');
% set(haGamma,'Children',[a(2:end);a(1)])

if lmeStatsBeta.pValue < 0.05
    hsigBeta = plot(haBeta,1.5,max(get(haBeta,'ylim')),'k*','markersize',10);
    legend([hsigBeta],['P<0.05'],'box','off');
end
if lmeStatsGamma.pValue < 0.05
    hsigGamma = plot(haGamma,1.5,max(get(haGamma,'ylim')),'k*','markersize',10);
    legend([hsigGamma],['P<0.05'],'box','off');
end
% save the figure to data location
saveas(hf,[pathname filename(1:end-4)])


% set(gco,'YData',get(gco,'YData')+22.5)
%% plot Bandpower
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/BP/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap);

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlightcorr = p_masked;

% the significant points could be marked with a red stars
highlight =double(highlight);
highlightcorr =double(highlightcorr);
highlight(highlight==0) = nan;
highlightcorr(highlightcorr==0) = nan;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
% for ic = 1:numel(Ls.PixelIdxList)
%     yy = Ls.PixelIdxList{ic};
%     if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%         highlight(Ls.PixelIdxList{ic}) = 0;
%     end
% end
% highlight =double(highlight);
% highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;


hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.time,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'k-','color',[0.5 0.5 0.5]);
hsigcorr = plot(Para.time,(max(y2plot(:))+0.1*range(y2plot(:)))*highlightcorr,'-','color',[0 0 0]);

plot([0,0],get(gca,'ylim'),'k--')
legend([hM.mainLine,hS.mainLine,hsig,hsigcorr], ...
    ['Intact'],['Scrambled'],['P<0.05'],['P<0.05 (corrected)']);
title(filename(1:end-4))
xlabel('Time relative to camera change (sec)')
ylabel('Normalised Power (a.u.)')
xlim([-0.5 1])
% ylim([min(y2plot(:))-0.3*range(y2plot(:)),max(y2plot(:))+0.4*range(y2plot(:))]);


% save the figure to data location
saveas(hf,[pathname filename(1:end-4)])

%% plot PLV
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/PLV/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap);

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.time,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'k*');
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.01 (corrected)']);
title(filename(1:end-4))
ylabel('Normalised Power (a.u)')
ylim([min(y2plot(:))-0.3*range(y2plot(:)),max(y2plot(:))+0.4*range(y2plot(:))]);


%% plot coherence

clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/COH/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap);

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.freq, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.freq, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.freq,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'k*');
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);
title(filename(1:end-4))
xlabel('Frequency (Hz)')
ylabel('Coherence (a.u)')
ylim([min(y2plot(:))-0.3*range(y2plot(:)),max(y2plot(:))+0.4*range(y2plot(:))]);


%% plot coherence coordinates correlation

clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/COH/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])

foi  = Para.freq>20 & Para.freq<30;

metricM = mean(allMetricM(:,foi,14),2);
metricS = mean(allMetricS(:,foi,14),2);

y = metricM;
x = allCoordinates(:,3);
p=polyfit(x,y,1);
yfit=polyval(p,x);

[R,P] = corrcoef(x,y);
[RS,PS] = corrcoef(allCoordinates(:,3),metricS);

hf1 = figure;
hp1 = plot(y,x,'*'); hold on;
hp2 = plot(yfit,x,'k-');

xlabel('Coherence')
ylabel('Z coordinates (mm)')

legend(hp2,['r = ' num2str(R(1,2)) ', p= ' num2str(P(1,2))])

rs = [];
ps = [];
for itime = 1:size(allMetricM,3)
    metricM = mean(allMetricM(:,foi,itime),2);
    metricS = mean(allMetricS(:,foi,itime),2);
    y = metricM;
    x = allCoordinates(:,3);
    [R,P] = corrcoef(x,y);
    rs(itime) = R(1,2);
    ps(itime) = P(1,2);
end
highlight =ones(size(ps));
highlight(ps>=0.05) = nan;

hf2 = figure;
plot(Para.timePT,rs);hold on;
hsig = plot(Para.timePT,max(rs)+0.3*range(rs)*highlight,'k*');

xlim([-0.5 1])
xlabel('Time relative to camera change')
ylabel('Correlation coefficient')

% save the figure to data location
saveas(hf1,[pathname filename(1:end-4)])
%% plot multivariate granger causality

clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/Granger/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])
% y2plot = yraw;
% se2plot = seraw;

% tMap(tMap==0)=nan;
% pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap);

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.freq, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.freq, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.freq,(max(y2plot(:))+0.1*range(y2plot(:)))*highlight,'k*');
plot([0,120],[0 0],'k--')
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);
title(filename(1:end-4))
ti = find(filename=='_');
ylabel([filename(1:ti-1) ' to Granger Index from ' filename(1:ti-1)])

xlim([1 120])
ylim([-1.5*max(abs(y2plot(:))),1.5*max(abs(y2plot(:)))])

% save the figure to data location
saveas(hf,[pathname filename(1:end-4)])

%% plot PSI
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/PSI/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])
y2plot = yraw;
tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap);

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.time,(max(y2plot(:))+0.1*range(y2plot(:)))*highlight,'k*');
title(filename(1:end-4))
xlim([-0.5 1])
ylim([-1.5*max(abs(y2plot(:))),1.5*max(abs(y2plot(:)))])
ti = find(filename=='_');
xlabel('Time relative to camera change (sec)')
ylabel([filename(1:ti-1) ' to  PSI  from ' filename(1:ti-1)])
plot([0,0],get(gca,'ylim'),'k--')

legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);

% save the figure to data location
saveas(hf,[pathname filename(1:end-4)])
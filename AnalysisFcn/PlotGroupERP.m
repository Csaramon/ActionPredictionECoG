%% section1: plot ERP
clear;

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

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
% ylim([ -3 3])
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')

legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);

% save the figure to data location
saveas(hf,[pathname filename(1:end-4)])

%% section2: plot Power Spectrum
clear;

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

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

% generete the main figure and specify the displaying style
hf = figure('Units','centimeters','Position', [8 6 16 12]);
set(gca,'linewidth',3,'FontSize',16)
hold on
xlim([0 120])


yyaxis left
set(gca,'ycolor','k')
ylabel('Low Frequency Power (a.u.)')
hlowM = shadedErrorBar(Para.freq(lowFind), y2plot(1,lowFind),se2plot(1,lowFind),{'color',[255 106 106]/255,'linewidth',1,'linestyle','-'},1);
hlowS = shadedErrorBar(Para.freq(lowFind), y2plot(2,lowFind), se2plot(2,lowFind),{'color',[30 144 255]/255,'linewidth',1,'linestyle','-'},1);
hlowS.mainLine.Marker = 'none';

hlowsig = plot(Para.freq(lowFind),(max(y2plot(:,lowFind),[],'all')+0.2*range(y2plot(:,lowFind),'all'))*highlight(lowFind),'-','color',[0.5 0.5 0.5],'LineWidth',3);
hlowsigcorr = plot(Para.freq(lowFind),(max(y2plot(:,lowFind),[],'all')+0.1*range(y2plot(:,lowFind),'all'))*highlightcorr(lowFind),'-','color',[0 0 0],'LineWidth',3);

yyaxis right
set(gca,'ycolor','k')
ylabel('High Frequency Power (a.u.)')
hhighM = shadedErrorBar(Para.freq(highFind), y2plot(1,highFind),se2plot(1,highFind),{'color',[255 106 106]/255,'linewidth',1,'linestyle','-'},1);
hhighS = shadedErrorBar(Para.freq(highFind), y2plot(2,highFind), se2plot(2,highFind),{'color',[30 144 255]/255,'linewidth',1,'linestyle','-'},1);
hhighS.mainLine.Marker = 'none';

hhighsig = plot(Para.freq(highFind),(max(y2plot(:,highFind),[],'all')+0.2*range(y2plot(:,highFind),'all'))*highlight(highFind),'-','color',[0.5 0.5 0.5],'LineWidth',3);
hhighsigcorr = plot(Para.freq(highFind),(max(y2plot(:,highFind),[],'all')+0.1*range(y2plot(:,highFind),'all'))*highlightcorr(highFind),'-','color',[0 0 0],'LineWidth',3);


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
    'Xlim',[0.5 2.5],'XTick', [1,2],'XTickLabel',{'I','S'},...
    'YTick',[],'FontSize',16,'NextPlot','add');
ylabel(haBeta,'β Power')
haGamma = axes(hf,'Position',[0.7 0.5 0.15 0.25],...
    'Xlim',[0.5 2.5],'XTick', [1,2],'XTickLabel',{'I','S'},...
    'YTick',[],'FontSize',16,'NextPlot','add');
ylabel(haGamma,'γ Power')


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
    plot(haBeta,[1,2],[ybetaM,ybetaS],'-','Color',[subColor(double(isub),:) 0.25],'LineWidth',3)
    
    ygammaM = mean(YcorrGamma(lmeTBL.Cond=='1' & lmeTBL.Sub==isub));
    ygammaS = mean(YcorrGamma(lmeTBL.Cond=='2' & lmeTBL.Sub==isub));
    plot(haGamma,[1,2],[ygammaM,ygammaS],'-','Color',[subColor(double(isub),:) 0.25],'LineWidth',3)
    
    ysubbetaM = [ysubbetaM,ybetaM];
    ysubbetaS = [ysubbetaS,ybetaS];
    ysubgammaM = [ysubgammaM,ygammaM];
    ysubgammaS = [ysubgammaS,ygammaS];
end

% plot mean value
plot(haBeta,[1,2],[nanmean(ysubbetaM),nanmean(ysubbetaS)],'-','Color',[0.4 0.4 0.4 1],'LineWidth',3)
set(haBeta,'LineWidth',3)
% a = get(haBeta,'Children');
% set(haBeta,'Children',[a(2:end);a(1)])
plot(haGamma,[1,2],[nanmean(ysubgammaM),nanmean(ysubgammaS)],'-','Color',[0.4 0.4 0.4 1],'LineWidth',3)
set(haGamma,'LineWidth',3)
% a = get(haGamma,'Children');
% set(haGamma,'Children',[a(2:end);a(1)])

if lmeStatsBeta.pValue < 0.05
    hsigBeta = plot(haBeta,1.5,max(get(haBeta,'ylim')),'k*','markersize',8);
    legend([hsigBeta],['P<0.05'],'box','off');
end
if lmeStatsGamma.pValue < 0.05
    hsigGamma = plot(haGamma,1.5,max(get(haGamma,'ylim')),'k*','markersize',8);
    legend([hsigGamma],['P<0.05'],'box','off');
end
% save the figure to data location
hf.Renderer = 'painters';
printeps(hf,[pathname filename(1:end-4)])


% set(gco,'YData',get(gco,'YData')+22.5)

%% section2-2: plot Bandpower
clear;

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

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

% generete the main figure and specify the displaying style
hf = figure('Units','centimeters','Position', [8 6 8 6]);
set(gca,'linewidth',1.5,'Units','centimeters','position',[1.5 1 6 4.5])
hold on

hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255,'linewidth',1.5,'linestyle','-'},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255,'linewidth',1.5,'linestyle','-'},1);
hsig = plot(Para.time,(max(y2plot(:))+0.2*range(y2plot(:)))*highlight,'k-','color',[0.5 0.5 0.5],'linewidth',1.5);
hsigcorr = plot(Para.time,(max(y2plot(:))+0.1*range(y2plot(:)))*highlightcorr,'-','color',[0 0 0],'linewidth',1.5);

ylim(get(gca,'ylim'))
plot([0,0],get(gca,'ylim'),'k--','linewidth',1)
legend([hM.mainLine,hS.mainLine,hsig,hsigcorr], ...
    ['Intact'],['Scrambled'],['P<0.05'],['P<0.05 (corrected)'] ...
    ,'box','off','NumColumns',2, ...
    'Position',[0.4 0.6 0.5 0.2]);

title(filename(1:end-4))
xlabel('Time relative to camera change (sec)')
ylabel('Normalised Power (a.u.)')
xlim([-0.5 1])

% save the figure to data location
hf.Renderer = 'painters';
set(findall(hf,'-property','FontSize'),'FontSize',10)
set(findall(hf,'-property','FontName'),'FontName','Arial')
set(findall(hf,'-property','FontWeight'),'FontWeight','Normal')
printeps(hf,[pathname filename(1:end-4)])

%% section3: plot PLV
clear;

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

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


%% section4: plot coherence
clear;

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

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

% generete the main figure and specify the displaying style
hf = figure('Units','centimeters','Position', [8 6 16 12]);
set(gca,'linewidth',3,'FontSize',20)

hold on
hM = shadedErrorBar(Para.freq, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255,'linewidth',3,'linestyle','-'},1);
hS = shadedErrorBar(Para.freq, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255,'linewidth',3,'linestyle','-'},1);
hsig = plot(Para.freq,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'-','color',[0 0 0],'linewidth',3);
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)'] ...
        ,'box','off','NumColumns',2);
title(filename(1:end-4))
xlabel('Frequency (Hz)')
ylabel('Coherence (a.u)')
ylim([min(y2plot(:))-0.3*range(y2plot(:)),max(y2plot(:))+0.4*range(y2plot(:))]);


% save the figure to data location
hf.Renderer = 'painters';
printeps(hf,[pathname filename(1:end-4)])


%% section4-2: plot freq coherence

clear ;

% mid time point of choosen time window
TimePoint = 0.5;

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

if ~filterindex
    return
end

% load in TF Coherence data
load([pathname filename])

% Choose time point
toi = Para.timePT > TimePoint-0.05 & Para.timePT < TimePoint+0.05;
MetricM = squeeze(nanmean(allMetricM(:,:,toi),3));
MetricS = squeeze(nanmean(allMetricS(:,:,toi),3));
x2plot = Para.freq;

% initialize result variable
tMap = nan(1,size(MetricM,2));
pMap = nan(1,size(MetricM,2));
y2plot = nan(2,size(MetricM,2));
se2plot = nan(2,size(MetricM,2));
yraw = nan(2,size(MetricM,2));
seraw = nan(2,size(MetricM,2));

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

y2plot = yraw;
se2plot = seraw;

% FDR  correction (only for 2 to 90 Hz)
corrFreqInd = Para.freq>2 & Para.freq< 90;
[p_fdr, p_masked] = fdr(pMap(corrFreqInd), 0.05,'Parametric');
highlight = pMap<=p_fdr;

highlight =double(highlight);
highlight(highlight==0) = nan;

% generete the main figure and specify the displaying style
hf = figure('Units','centimeters','Position', [8 6 16 12]);
set(gca,'linewidth',3,'FontSize',20)

hold on
hM = shadedErrorBar(x2plot, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255,'linewidth',3,'linestyle','-'},1);
hS = shadedErrorBar(x2plot, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255,'linewidth',3,'linestyle','-'},1);

hsig = plot(x2plot,(max(y2plot(:))+0.2*range(y2plot(:)))*highlight,'k-','linewidth',3);
% ignore frequencies around 50Hz line noise
rectangle('Position',[45 min(get(gca,'ylim'))+0.1*range(get(gca,'ylim')) ...
    10 0.8*range(get(gca,'ylim'))], ...
    'FaceColor','w','EdgeColor','w')

xlim([0 90])
plot([0,0],get(gca,'ylim'),'k--','linewidth',1.5);
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['p<0.05 (corrected)'] ...
    ,'box','off','NumColumns',2);

ti = find(filename=='_');
title([filename(1:ti-1) '--' filename(ti+1:end-4) '(' ...
    num2str(TimePoint-0.5*Para.timeWin) '-' ...
    num2str(TimePoint+0.5*Para.timeWin) 's)'])
xlabel('Frequency (Hz)')
ylabel('Imaginary Coherence')

% save the figure to data location
hf.Renderer = 'painters';
printeps(hf,[pathname filename(1:end-4) num2str(TimePoint) 's'])


%% section4-3: plot time coherence

clear ;

% mid time point of choosen time window
% FreqRange = [20 30];
FreqRange = [60 90];

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

if ~filterindex
    return
end

% load in TF Coherence data
load([pathname filename])

% Choose time point
foi  = Para.freq>min(FreqRange) & Para.freq<max(FreqRange);
MetricM = squeeze(mean(allMetricM(:,foi,:),2));
MetricS = squeeze(mean(allMetricS(:,foi,:),2));
x2plot = Para.timePT;

% initialize result variable
tMap = nan(1,size(MetricM,2));
pMap = nan(1,size(MetricM,2));
y2plot = nan(2,size(MetricM,2));
se2plot = nan(2,size(MetricM,2));
yraw = nan(2,size(MetricM,2));
seraw = nan(2,size(MetricM,2));

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

y2plot = yraw;
se2plot = seraw;

% FDR  correction (only for -0.5 to 1s relative to camera change)
corrTimeInd = Para.timePT>-0.55 & Para.timePT< 1.05;
[p_fdr, p_masked] = fdr(pMap(corrTimeInd), 0.05,'Parametric');
highlight = pMap<=p_fdr;

highlight =double(highlight);
highlight(highlight==0) = nan;

% generete the main figure and specify the displaying style
hf = figure('Units','centimeters','Position', [8 6 16 12]);
set(gca,'linewidth',3,'FontSize',20)

hold on
hM = shadedErrorBar(x2plot, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255,'linewidth',3,'linestyle','-'},1);
hS = shadedErrorBar(x2plot, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255,'linewidth',3,'linestyle','-'},1);

hsig = plot(x2plot,(max(y2plot(:))+0.2*range(y2plot(:)))*highlight,'k-','linewidth',3);

xlim([-0.5 1])
plot([0,0],get(gca,'ylim'),'k--','linewidth',1.5);
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['p<0.05 (corrected)'] ...
    ,'box','off');

ti = find(filename=='_');
title([filename(1:ti-1) '--' filename(ti+1:end-4) '(' ...
    num2str(min(FreqRange)) '-' ...
    num2str(max(FreqRange)) 'Hz)'])
xlabel('Time relative to camera change (sec)')
ylabel('Imaginary Coherence')

% save the figure to data location
hf.Renderer = 'painters';
printeps(hf,[pathname filename(1:end-4) ...
    num2str(min(FreqRange)) '_' num2str(max(FreqRange)) 'Hz'])


%% section5: plot freq  non-parametric granger causality

clear ;

% mid time point of choosen time window
TimePoint = 0.5;

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

if ~filterindex
    return
end

% load in TF Coherence data
load([pathname filename])

% Choose time point
toi = Para.timePT > TimePoint-0.05 & Para.timePT < TimePoint+0.05;
MetricM = squeeze(nanmean(allMetricM(:,:,toi),3));
MetricS = squeeze(nanmean(allMetricS(:,:,toi),3));
x2plot = Para.freq;

% initialize result variable
tMap = nan(1,size(MetricM,2));
pMap = nan(1,size(MetricM,2));
y2plot = nan(2,size(MetricM,2));
se2plot = nan(2,size(MetricM,2));
yraw = nan(2,size(MetricM,2));
seraw = nan(2,size(MetricM,2));

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

y2plot = yraw;
se2plot = seraw;

% FDR  correction (only for 2 to 90 Hz)
corrFreqInd = Para.freq>2 & Para.freq< 90;
[p_fdr, p_masked] = fdr(pMap(corrFreqInd), 0.05,'Parametric');
highlight = pMap<=p_fdr;

highlight =double(highlight);
highlight(highlight==0) = nan;

% generete the main figure and specify the displaying style
hf = figure('Units','centimeters','Position', [8 6 12 9]);
set(gca,'linewidth',1.5,'Units','centimeters','position',[2 1.5 6 4.5])

hold on
hM = shadedErrorBar(x2plot, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255,'linewidth',1.5,'linestyle','-'},1);
hS = shadedErrorBar(x2plot, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255,'linewidth',1.5,'linestyle','-'},1);

hsig = plot(x2plot,(max(y2plot(:))+0.2*range(y2plot(:)))*highlight,'k-','linewidth',1.5);
% ignore frequencies around 50Hz line noise
rectangle('Position',[45 min(get(gca,'ylim'))+0.1*range(get(gca,'ylim')) ...
    10 0.8*range(get(gca,'ylim'))], ...
    'FaceColor','w','EdgeColor','w')

xlim([0 90])
plot([0,0],get(gca,'ylim'),'k--','linewidth',1);
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['p<0.05 (corrected)'] ...
    ,'box','off','NumColumns',2, ...
    'Units','centimeters', 'Position',[2.5 4 4 1.5]);

ti = find(filename=='_');
title([filename(1:ti-1) '--' filename(ti+1:end-4) '(' ...
    num2str(TimePoint-0.5*Para.timeWin) '-' ...
    num2str(TimePoint+0.5*Para.timeWin) 's)'])
xlabel('Frequency (Hz)')
ylabel({'GC',['from ' filename(1:ti-1) ' to ' filename(ti+1:end-4)]})

% save the figure to data location
hf.Renderer = 'painters';
set(findall(hf,'-property','FontSize'),'FontSize',10)
set(findall(hf,'-property','FontName'),'FontName','Arial')
set(findall(hf,'-property','FontWeight'),'FontWeight','Normal')
printeps(hf,[pathname filename(1:end-4) num2str(TimePoint) 's'])


%% section6: plot time PSI

clear ;

% mid time point of choosen time window
% FreqRange = [20 30];
FreqRange = [60 90];

T = which('PlotGroupERP');
[runpath,~,~] = fileparts(T);
[filename, pathname, filterindex] = uigetfile([runpath(1:end-11) '/Results/*.mat']);

if ~filterindex
    return
end

% load in TF Coherence data
load([pathname filename])

% Choose time point
foi  = Para.freq>min(FreqRange) & Para.freq<max(FreqRange);
MetricM = squeeze(mean(allMetricM(:,foi,:),2));
MetricS = squeeze(mean(allMetricS(:,foi,:),2));
x2plot = Para.timePT;

% initialize result variable
tMap = nan(1,size(MetricM,2));
pMap = nan(1,size(MetricM,2));
y2plot = nan(2,size(MetricM,2));
se2plot = nan(2,size(MetricM,2));
yraw = nan(2,size(MetricM,2));
seraw = nan(2,size(MetricM,2));

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

y2plot = yraw;
se2plot = seraw;

% FDR  correction (only for -0.5 to 1s relative to camera change)
corrTimeInd = Para.timePT>-0.55 & Para.timePT< 1.05;
[p_fdr, p_masked] = fdr(pMap(corrTimeInd), 0.05,'Parametric');
highlight = pMap<=p_fdr;

highlight =double(highlight);
highlight(highlight==0) = nan;

% generete the main figure and specify the displaying style
hf = figure('Units','centimeters','Position', [8 6 8 6]);
set(gca,'linewidth',1.5,'Units','centimeters','position',[1.5 1 6 4.5])
hold on

hM = shadedErrorBar(x2plot, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255,'linewidth',1.5,'linestyle','-'},1);
hS = shadedErrorBar(x2plot, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255,'linewidth',1.5,'linestyle','-'},1);

hsig = plot(x2plot,(max(y2plot(:))+0.2*range(y2plot(:)))*highlight,'k-','linewidth',1.5);

xlim([-0.5 1])
plot([0,0],get(gca,'ylim'),'k--','linewidth',1);
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['p<0.05 (corrected)'] ...
    ,'box','off','NumColumns',2, ...
    'Position',[0.4 0.6 0.5 0.2]);

ti = find(filename=='_');
title([filename(1:ti-1) '--' filename(ti+1:end-4) '(' ...
    num2str(min(FreqRange)) '-' ...
    num2str(max(FreqRange)) 'Hz)'])
xlabel('Time relative to camera change (sec)')
ylabel({'PSI',['←' filename(1:ti-1) ' from ' filename(ti+1:end-4) '→']})

% save the figure to data location
hf.Renderer = 'painters';
set(findall(hf,'-property','FontSize'),'FontSize',10)
set(findall(hf,'-property','FontName'),'FontName','Arial')
set(findall(hf,'-property','FontWeight'),'FontWeight','Normal')
printeps(hf,[pathname filename(1:end-4) ...
    num2str(min(FreqRange)) '_' num2str(max(FreqRange)) 'Hz'])


%% plot coherence coordinates correlation

clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Documents/ActionPredictionECoG/Results/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])

foi  = Para.freq>20 & Para.freq<30;

metricM = mean(allMetricM(:,foi,14),2);
metricS = mean(allMetricS(:,foi,14),2);

% linear fit
y = metricM;
x = allCoordinates(:,3);
p=polyfit(x,y,1);
yfit=polyval(p,x);

ys = metricS;
x = allCoordinates(:,3);
ps=polyfit(x,ys,1);
yfits=polyval(ps,x);

% plot correlation 
[R,P] = corrcoef(x,y);
[Rs,Ps] = corrcoef(x,ys);
R = round(R*100)/100;P = round(P*100)/100;
Rs = round(Rs*100)/100;Ps = round(Ps*100)/100;

% generete the main figure and specify the displaying style
hf1 = figure('Units','centimeters','Position', [8 6 8 6]);
set(gca,'linewidth',1.5,'Units','centimeters','position',[1.5 1 6 4.5])
hold on

hp1 = plot(y,x,'.','color',[255 106 106]/255,'markersize',10);
hp1s = plot(ys,x,'.','color',[30 144 255]/255,'markersize',10);
hp2 = plot(yfit,x,'-','color',[255 106 106]/255,'linewidth',1);
hp2s = plot(yfits,x,'-','color',[30 144 255]/255,'linewidth',1);

ylim([min(x)-5 max(x)+5])
xlabel('Beta Coherence')
ylabel('Z coordinates (mm)')

legend([hp2,hp2s],['Intact: r = ' num2str(R(1,2)) ', p= ' num2str(P(1,2))],...
    ['Scrambled: r = ' num2str(Rs(1,2)) ', p= ' num2str(Ps(1,2))],...
    'box','off','Position',[0.4 0.6 0.5 0.2]);

% plot time variant correlation
for itime = 1:size(allMetricM,3)
    metricM = mean(allMetricM(:,foi,itime),2);
    metricS = mean(allMetricS(:,foi,itime),2);
    y = metricM;
    ys = metricS;
    x = allCoordinates(:,3);
    [R,P] = corrcoef(x,y);
    rts(itime) = R(1,2);
    pts(itime) = P(1,2);
    
    [Rs,Ps] = corrcoef(x,ys);
    rts2(itime) = Rs(1,2);
    pts2(itime) = Ps(1,2);
end
highlight =ones(size(pts));
highlight(pts>=0.05) = nan;
highlight2 =ones(size(pts2));
highlight2(pts2>=0.05) = nan;

% generete the main figure and specify the displaying style
hf2 = figure('Units','centimeters','Position', [8 6 8 6]);
set(gca,'linewidth',1.5,'Units','centimeters','position',[1.5 1 6 4.5])
hold on

hl1 = plot(Para.timePT,rts,'color',[255 106 106]/255,'linewidth',1.5);
hl1s = plot(Para.timePT,rts2,'color',[30 144 255]/255,'linewidth',1.5);
hsig = plot(Para.timePT,max(rts)+0.3*range(rts)*highlight,'k*','linewidth',1.5);
hsigs = plot(Para.timePT,max(rts2)+0.3*range(rts2)*highlight2,'k*','linewidth',1.5);

legend([hl1,hl1s],['Intact'],['Scrambled'], ...
    'box','off','Position',[0.4 0.6 0.5 0.2]);

xlim([-0.5 1])
xlabel('Time relative to camera change')
ylabel('Correlation coefficient')


% save the figure to data location
hf1.Renderer = 'painters';
set(findall(hf1,'-property','FontSize'),'FontSize',10)
set(findall(hf1,'-property','FontName'),'FontName','Arial')
set(findall(hf1,'-property','FontWeight'),'FontWeight','Normal')
printeps(hf1,[pathname filename(1:end-4)])
hf2.Renderer = 'painters';
set(findall(hf2,'-property','FontSize'),'FontSize',10)
set(findall(hf2,'-property','FontName'),'FontName','Arial')
set(findall(hf2,'-property','FontWeight'),'FontWeight','Normal')
printeps(hf2,[pathname filename(1:end-4) 'temporal'])
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

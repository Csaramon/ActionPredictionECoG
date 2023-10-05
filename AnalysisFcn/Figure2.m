%% User defined data path
dataPath = 'C:\Users\qin2\Documents\ActionPredictionECoG\Data\ECoGData\';
addpath('C:\Users\qin2\Documents\MATLAB\toolbox\fieldtrip-20210418') % path to fieldtrip toolbox


%% Calculate PSD
allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub

% Region of Interest include:
ROIIndex = {[1,2],[63,64],[51,52]};
ROIAtlas = {'fsAnatomyMacro.nii','fsAnatomyMacro.nii','fsAnatomyMacro.nii'};
ROIText = {'Precentral','SupraMarginal','MiddleOccipitalGyrus'};

roiDist = 1; % maximum distance between electrodes and ROI voxels

for iatlas = [1,2,3]
    
    % load altlas infomation
    aparc_nii = load_nifti(['Atlas' filesep ROIAtlas{iatlas}]);
    aparc_vol = round(aparc_nii.vol);
    allroi_index = [];
    %%%------ For maxprobabilistic map ------%%%
    for iroi = ROIIndex{iatlas}
        temp_indx = find(aparc_vol==iroi);
        allroi_index = [allroi_index;temp_indx];
    end
    
    % atlas region coordinates
    [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
    ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
    ras_coordiantes = ras_coordiantes(1:3,:)';
    aparc_coordiantes = ras_coordiantes;
    
    % Initialize variables
    allMetricM = [];
    allMetricS = [];
    subIndexM = [];
    subIndexS = [];
    elecIndexM = [];
    elecIndexS =[];
    trlIndexM = [];
    trlIndexS = [];
    nelec = 1;
    Para.elecposMNI = [];
    Para.time = [];
    Para.freq = [];
    
    % -------- Subject Level -------- %
    for isub = 1:numel(allsub)
        subname = allsub{isub};
        fprintf(['\n Currently calculating subject: ' subname])
        
        %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist([dataPath subname 'LARER_rerefDataInROI.mat'],'file')
            a = load([dataPath subname 'LARER_rerefDataInROI.mat']);
        end
        c = fieldnames(a);
        rerefData = a.(c{1});
        
        % choose ROI electrodes according to MNI coordinates
        elecposMNI = rerefData.elec.elecposMNI;
        tempdev = pdist2(elecposMNI,aparc_coordiantes);
        [it,~] = find(tempdev <=roiDist);
        ROIelec = unique(it);
        
        % skip if no electrode pair
        if ~any(ROIelec)
            continue
        end
        
        if exist([dataPath subname '_eventdata.mat'],'file')
            a = load([dataPath subname '_eventdata']);
        end
        c = fieldnames(a);
        camInfo = a.(c{1});
        
        fs = rerefData.fsample;
        camInfo(:,4) = num2cell(cell2mat(camInfo(:,4))./1000);
        
        
        % time-frequency decomposition (multi-taper)
        cfg              = [];
        cfg.channel = ROIelec;
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.foi          = 32:2:120;
        cfg.t_ftimwin  = 0.5*ones(size(cfg.foi));%14./cfg.foi;
        cfg.tapsmofrq  = 6*ones(size(cfg.foi));%0.3 *cfg.foi; % tapers=2*tw*fw-1
        cfg.toi          = min(rerefData.time{1}):0.1:max(rerefData.time{1});
        %                                 cfg.precision = 'single';
        cfg.pad='nextpow2';
        
        
        ft_warning off
        freq2 = ft_freqanalysis(cfg,rerefData);
        
        cfg.foi          = 2:1:30;
        cfg.t_ftimwin  = 1*ones(size(cfg.foi));
        cfg.taper      =  'hanning';
        freq = ft_freqanalysis(cfg,rerefData);
        
        clear rerefData
        
        allTS = cat(2,freq.powspctrm,freq2.powspctrm);
        
        % predefine variable for the reconcantenate time series
        tsI = [];
        tsS = [];
        %%%%-------- reconcatenate data movie clip-wise and condition-wise --------%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for im = 1:size(camInfo,1)
            
            timeInd = freq.time >= camInfo{im,4} & freq.time <= (camInfo{im,4}+camInfo{im,6});
            % extract data of each movie
            tmpTS = allTS(:,:,timeInd);
            
            
            % concatenate data according to condition
            if camInfo{im,1}==0
                tsI=cat(3,tsI,mean(tmpTS,3,'omitnan'));
            elseif camInfo{im,1}==1
                tsS=cat(3,tsS,mean(tmpTS,3,'omitnan'));
            end
            
            
        end
        
        % generate index of Subject Electrode and Trial for LME
        metricM = mean(tsI,3);
        metricS = mean(tsS,3);
        
        elecIndexM =  cat(1,elecIndexM,[nelec:nelec+numel(ROIelec)-1]');
        subIndexM = cat(1,subIndexM,repmat(isub,numel(ROIelec),1));
        
        elecIndexS = cat(1,elecIndexS,[nelec:nelec+numel(ROIelec)-1]');
        subIndexS = cat(1,subIndexS,repmat(isub,numel(ROIelec),1));
        
        % concontenate all trial responses in chosen ROI
        allMetricM = cat(1,allMetricM,metricM);
        allMetricS = cat(1,allMetricS,metricS);
        
        
        Para.elecposMNI = [Para.elecposMNI;freq.elec.elecposMNI(ROIelec,:)];
        Para.freq = cat(2,freq.freq,freq2.freq);
        
        nelec = nelec+numel(ROIelec);
        
        
        
    end
    
    
    %% Fit power spectrum LME
    
    CondIndexM = ones(size(subIndexM));
    CondIndexS = 2*ones(size(subIndexS));
    
    tMap = zeros(1,size(allMetricM,2));
    pMap = zeros(1,size(allMetricM,2));
    y2plot = zeros(2,size(allMetricM,2));
    se2plot = zeros(2,size(allMetricM,2));
    yraw = zeros(2,size(allMetricM,2));
    seraw = zeros(2,size(allMetricM,2));
    
    strlen = 0;
    % TFR point wise LME
    
    for itime = 1:size(allMetricM,2)
        
        s = ['Calculating tf point: ' num2str(itime) '/' num2str(size(allMetricM,2)) 'times'];
        strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
        strlen = strlentmp - strlen;
        
        frameDataM = double(squeeze(allMetricM(:,itime)));
        % skip nan point
        if ~any(frameDataM,'all')
            continue
        end
        frameDataS = double(squeeze(allMetricS(:,itime)));
        
        lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
            [elecIndexM;elecIndexS], 'VariableNames',{'Y','Cond','Sub','Elec'});
        lmeTBL.Cond = nominal(lmeTBL.Cond);
        lmeTBL.Sub = nominal(lmeTBL.Sub);
        lmeTBL.Elec = nominal(lmeTBL.Elec);
        lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
        
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
    
    
    %% Plot power spectrum results
    
    tMap(tMap==0)=nan;
    pMap(tMap==0)=nan;
    
    lowFind = Para.freq <= 30;
    highFind = Para.freq > 30;
    % rescale the results to compensate for the 1/f trend of the PSD
    y2plot = y2plot.*Para.freq.^2;
    se2plot = se2plot.*Para.freq.^2;
    
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
    hf = figure('Units','centimeters','Position', [8 6 32 27]);
    set(gca,'linewidth',5,'FontSize',28,'FontName','Arial','InnerPosition',[0.11,0.11,0.78,0.82])
    hold on
    xlim([0 120])
    
    
    yyaxis left
    set(gca,'ycolor','k')
    ylabel('Low Frequency Power (a.u.)')
    hlowM = shadedErrorBar(Para.freq(lowFind), y2plot(1,lowFind),se2plot(1,lowFind),{'color',[255 106 106]/255,'linewidth',1.5,'linestyle','-'},1);
    hlowS = shadedErrorBar(Para.freq(lowFind), y2plot(2,lowFind), se2plot(2,lowFind),{'color',[30 144 255]/255,'linewidth',1.5,'linestyle','-'},1);
    hlowS.mainLine.Marker = 'none';
    
    hlowsig = plot(Para.freq(lowFind),(max(y2plot(:,lowFind),[],'all')+0.2*range(y2plot(:,lowFind),'all'))*highlight(lowFind),'-','color',[0.5 0.5 0.5],'LineWidth',5);
    hlowsigcorr = plot(Para.freq(lowFind),(max(y2plot(:,lowFind),[],'all')+0.1*range(y2plot(:,lowFind),'all'))*highlightcorr(lowFind),'-','color',[0 0 0],'LineWidth',5);
    
    yyaxis right
    set(gca,'ycolor','k')
    ylabel('High Frequency Power (a.u.)')
    hhighM = shadedErrorBar(Para.freq(highFind), y2plot(1,highFind),se2plot(1,highFind),{'color',[255 106 106]/255,'linewidth',1.5,'linestyle','-'},1);
    hhighS = shadedErrorBar(Para.freq(highFind), y2plot(2,highFind), se2plot(2,highFind),{'color',[30 144 255]/255,'linewidth',1.5,'linestyle','-'},1);
    hhighS.mainLine.Marker = 'none';
    
    hhighsig = plot(Para.freq(highFind),(max(y2plot(:,highFind),[],'all')+0.2*range(y2plot(:,highFind),'all'))*highlight(highFind),'-','color',[0.5 0.5 0.5],'LineWidth',5);
    hhighsigcorr = plot(Para.freq(highFind),(max(y2plot(:,highFind),[],'all')+0.1*range(y2plot(:,highFind),'all'))*highlightcorr(highFind),'-','color',[0 0 0],'LineWidth',5);
    
    
    title(ROIText{iatlas})
    xlabel('Frequency (Hz)')
    
    legend([hlowM.mainLine,hlowS.mainLine,hlowsig,hlowsigcorr], ...
        ['Intact'],['Scrambled'],['p<0.05'],['p<0.05 (corrected)'],'box','off','NumColumns',2);
    
    
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
        plot(haBeta,[1,2],[ybetaM,ybetaS],'-','Color',[subColor(double(isub),:) 0.25],'LineWidth',5)
        
        ygammaM = mean(YcorrGamma(lmeTBL.Cond=='1' & lmeTBL.Sub==isub));
        ygammaS = mean(YcorrGamma(lmeTBL.Cond=='2' & lmeTBL.Sub==isub));
        plot(haGamma,[1,2],[ygammaM,ygammaS],'-','Color',[subColor(double(isub),:) 0.25],'LineWidth',5)
        
        ysubbetaM = [ysubbetaM,ybetaM];
        ysubbetaS = [ysubbetaS,ybetaS];
        ysubgammaM = [ysubgammaM,ygammaM];
        ysubgammaS = [ysubgammaS,ygammaS];
    end
    
    % plot mean value
    plot(haBeta,[1,2],[nanmean(ysubbetaM),nanmean(ysubbetaS)],'-','Color',[0.4 0.4 0.4 1],'LineWidth',5)
    set(haBeta,'LineWidth',5)
    % a = get(haBeta,'Children');
    % set(haBeta,'Children',[a(2:end);a(1)])
    plot(haGamma,[1,2],[nanmean(ysubgammaM),nanmean(ysubgammaS)],'-','Color',[0.4 0.4 0.4 1],'LineWidth',5)
    set(haGamma,'LineWidth',5)
    % a = get(haGamma,'Children');
    % set(haGamma,'Children',[a(2:end);a(1)])
    
    if lmeStatsBeta.pValue < 0.05
        hsigBeta = plot(haBeta,1.5,max(get(haBeta,'ylim')),'k*','markersize',12,'linewidth',2);
        legend([hsigBeta],['p<0.05'],'box','off');
    end
    if lmeStatsGamma.pValue < 0.05
        hsigGamma = plot(haGamma,1.5,max(get(haGamma,'ylim')),'k*','markersize',12,'linewidth',2);
        legend([hsigGamma],['p<0.05'],'box','off');
    end
    
    % save the figure to data location
    % hf.Renderer = 'painters';
    % printeps(hf,[pathname filename(1:end-4)])
    
    % set(gco,'YData',get(gco,'YData')+22.5)
    
    
end
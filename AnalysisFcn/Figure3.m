%% User defined data path
dataPath = 'C:\Users\qin2\Documents\ActionPredictionECoG\Data\ECoGData\';
addpath('C:\Users\qin2\Documents\MATLAB\toolbox\fieldtrip-20210418') % path to fieldtrip toolbox


%% Calculate Band power
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
        
        if exist([dataPath subname '_eventdata.mat'],'file')
            a = load([dataPath subname '_eventdata']);
        end
        c = fieldnames(a);
        camInfo = a.(c{1});
        
        % choose ROI electrodes according to MNI coordinates
        elecposMNI = rerefData.elec.elecposMNI;
        tempdev = pdist2(elecposMNI,aparc_coordiantes);
        [it,~] = find(tempdev <=roiDist);
        ROIelec = unique(it);
        
        % Seperate Trials
        preTrigger = 1.3; % -0.5-0.75s  (2Hz 3cycle =0.75s half window)
        postTrigger = 1.8; % 1+0.75s (2Hz 3cycle =0.75s half window)
        
        numCams = size(camInfo,1);
        trl = nan(numCams,8);
        for cam = 1:numCams
            trl(cam,1:8) = [camInfo{cam,4}-(preTrigger*rerefData.fsample)...
                (camInfo{cam,4}+(postTrigger*rerefData.fsample)-1)...
                -preTrigger*rerefData.fsample...
                camInfo{cam,1}...
                camInfo{cam,3}...
                camInfo{cam,4}...
                camInfo{cam,5}...
                camInfo{cam,6}];
        end
        
        % construct fieldtrip trial data
        cfg = [];
        cfg.trl = trl;
        trlData = ft_redefinetrial(cfg,rerefData);
        
        timeWin = [-0.5 1];
        freqRange = [60 90] ; % [20 30] [60 90]
        
        % skip bad channels
        badChanInd = trlData.trial{1,1}(ROIelec,1)==0;
        ROIelec(badChanInd) = [];
        
        % skip if no electrode pair
        if ~any(ROIelec)
            continue
        end
        
        % seperate conditions
        cfg = [];
        cfg.trials = find(trlData.trialinfo(:,1)==0);
        trlDataM = ft_selectdata(cfg,trlData);
        
        cfg = [];
        cfg.trials = find(trlData.trialinfo(:,1)==1);
        trlDataS = ft_selectdata(cfg,trlData);
        
        clear trlData
        
        % multi taper
        cfg              = [];
        cfg.channel = ROIelec;
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.foi          = min(freqRange):max(freqRange);
        cfg.t_ftimwin  = 0.5*ones(size(cfg.foi));%14./cfg.foi;
        cfg.tapsmofrq  = 6*ones(size(cfg.foi));%0.3 *cfg.foi; % tapers=2*tw*fw-1
        cfg.toi          = min(timeWin):0.01:max(timeWin); % 'all';
        cfg.keeptrials = 'yes';
        cfg.pad='nextpow2';
        
        ft_warning off
        freqM = ft_freqanalysis(cfg,trlDataM);
        
        freqS = ft_freqanalysis(cfg,trlDataS);
        
        % raw power
        
        % normalise power
        tmpPowerM = freqM.powspctrm;
        tmpPowerS =freqS.powspctrm;

        metricM = tmpPowerM;
        metricM = nanmean(metricM,3);
        
        metricS = tmpPowerS;
        metricS = nanmean(metricS,3);
        
        
        %---------- elec level ----------%
        metricM = squeeze(mean(metricM,1));
        metricS = squeeze(mean(metricS,1));
        if size(metricM,1) > size(metricM,2)
            metricM = metricM';
            metricS =metricS';
        end
        
        elecIndexM =  cat(1,elecIndexM,[nelec:nelec+numel(ROIelec)-1]');
        subIndexM = cat(1,subIndexM,repmat(isub,numel(ROIelec),1));
        
        elecIndexS = cat(1,elecIndexS,[nelec:nelec+numel(ROIelec)-1]');
        subIndexS = cat(1,subIndexS,repmat(isub,numel(ROIelec),1));
        
        % concontenate all trial responses in chosen ROI
        allMetricM = cat(1,allMetricM,metricM);
        allMetricS = cat(1,allMetricS,metricS);
        
        
        Para.elecposMNI = [Para.elecposMNI;freqM.elec.elecposMNI(ROIelec,:)];
        Para.time = freqM.time;
        Para.freq = freqM.freq;
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

title(ROIText{iatlas})
xlabel('Time relative to camera change (sec)')
ylabel('Normalised Power (a.u.)')
xlim([-0.5 1])

end
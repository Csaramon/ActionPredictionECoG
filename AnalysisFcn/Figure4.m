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


seedIndex = [1 2 3];
searchIndex = [1 2 3];

allPair = nchoosek(seedIndex,2);
for iseed = seedIndex
    for isearch = searchIndex
        
        % skip redundant pairs
        if ~ismember([iseed,isearch],allPair,'rows')
            continue
        end
        
        % load altlas infomation
        aparc_nii = load_nifti(['Atlas' filesep ROIAtlas{iseed}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        for iroi = ROIIndex{iseed}
            temp_indx = find(aparc_vol==iroi);
            allroi_index = [allroi_index;temp_indx];
        end
        
        % atlas seed region coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        seed_coordiantes = ras_coordiantes;
        
        % load altlas infomation
        aparc_nii = load_nifti(['Atlas' filesep ROIAtlas{isearch}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        for iroi = ROIIndex{isearch}
            temp_indx = find(aparc_vol==iroi);
            allroi_index = [allroi_index;temp_indx];
        end
        
        % atlas search region coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        search_coordiantes = ras_coordiantes;
        
        
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
            
            p=0.05; % threshold for IVC
            timeWin = 1; % unit in second
            timeStep = 0.05; % unit in second
            
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
            
            cfg = [];
            cfg.demean = 'yes';
            cfg.detrend = 'yes';
            trlData = ft_preprocessing(cfg,trlData);
            
            % choose seed electrodes according to MNI coordinates
            elecposMNI = trlData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,seed_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            seedElec = unique(it);
            
            % choose search electrodes according to MNI coordinates
            tempdev = pdist2(elecposMNI,search_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            searchElec = unique(it);
            
            % skip bad channels
            badChanInd = trlData.trial{1,1}(seedElec,1)==0;
            seedElec(badChanInd) = [];
            
            badChanInd = trlData.trial{1,1}(searchElec,1)==0;
            searchElec(badChanInd) = [];
            
            % skip if no electrode pair
            if ~any(seedElec) | ~any(searchElec) | isempty(setdiff(seedElec,searchElec))
                continue
            end
            
            % seperate conditions
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==0);
            trlDataM = ft_selectdata(cfg,trlData);
            
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==1);
            trlDataS = ft_selectdata(cfg,trlData);
            
            timeRange = [min(trlData.time{1}) max(trlData.time{1})];
            
            % calculate pairwise COH
            Npair = 1;
            allChanCmb = [];
            duplicateInd = intersect(searchElec,seedElec);
            if ~isempty(duplicateInd)
                if numel(searchElec) > numel(seedElec)
                    searchElec = setdiff(searchElec,duplicateInd);
                else
                    seedElec = setdiff(seedElec,duplicateInd);
                end
            end
            
            in = 1;
            allcohM = [];
            allcohS = [];
            for itw = min(timeRange):timeStep:max(timeRange)-timeWin
                
                % choose time window
                cfg = [];
                cfg.latency = [itw itw+timeWin];
                trlDataMtmp = ft_selectdata(cfg,trlDataM);
                
                cfg = [];
                cfg.latency = [itw itw+timeWin];
                trlDataStmp = ft_selectdata(cfg,trlDataS);
                
                % calculate fourier spectrum using multi taper
                cfg            = [];
                cfg.output     = 'fourier';
                cfg.method     = 'mtmfft';
                %                                                             cfg.foilim     = [0 30];
                %                     cfg.foi          = 2:1:30;
                %                                         cfg.taper      =  'hanning';
                cfg.tapsmofrq  = 6;
                cfg.keeptrials = 'yes';
                
                freqM    = ft_freqanalysis(cfg, trlDataMtmp);
                
                freqS    = ft_freqanalysis(cfg, trlDataStmp);
                
                % calculate COH in Intact Condition
                cfg            = [];
                cfg.method     = 'coh';
                cfg.complex = 'absimag';
                cfg.channelcmb = {trlDataM.label(seedElec) trlDataM.label(searchElec)};
                COHM             = ft_connectivityanalysis(cfg, freqM);
                allcohM(:,:,in) = COHM.cohspctrm;
                % calculate COH in Scambled Condition
                cfg.channelcmb = {trlDataS.label(seedElec) trlDataS.label(searchElec)};
                COHS             = ft_connectivityanalysis(cfg, freqS);
                allcohS(:,:,in) = COHS.cohspctrm;
                
                timePt(in) = itw+0.5*timeWin;
                in = in+1;
            end
            
            allChanCmb = [allChanCmb;COHM.labelcmb];
            
            % count for pairs of eletrodes
            Npair = Npair + size(COHM.labelcmb,1);
            
            % generate index for Subject Electrode and Trial
            metricM = allcohM;
            metricS = allcohS;
            
            elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
            subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
            
            elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
            subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,metricM);
            allMetricS = cat(1,allMetricS,metricS);
            
            
            Para.chanCMB{isub} = allChanCmb;
            Para.timeWin = timeWin;
            Para.timeStep = timeStep;
            Para.timePT = timePt;
            Para.freq = freqM.freq;
            nelec = nelec+Npair-1;
            
        end
        
        
        %% section2-2: calculate Coherence (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CondIndexM = ones(size(subIndexM));
        CondIndexS = 2*ones(size(subIndexS));
        
        tMap = zeros(size(allMetricS,2),size(allMetricS,3));
        pMap = zeros(size(allMetricS,2),size(allMetricS,3));
        y2plot = zeros(2,size(allMetricS,2),size(allMetricS,3));
        se2plot = zeros(2,size(allMetricS,2),size(allMetricS,3));
        yraw = zeros(2,size(allMetricS,2),size(allMetricS,3));
        seraw = zeros(2,size(allMetricS,2),size(allMetricS,3));
        
        strlen = 0;
        % TFR point wise LME
        for ifreq = 1:size(allMetricM,2)
            for itime = 1:size(allMetricM,3)
                
                s = ['Calculating tf point: ' num2str(ifreq) '/' num2str(size(allMetricM,2)) 'freqs, ' ...
                    num2str(itime) '/' num2str(size(allMetricM,3)) 'times'];
                strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
                strlen = strlentmp - strlen;
                
                frameDataM = double(squeeze(allMetricM(:,ifreq,itime)));
                % skip nan point
                if ~any(frameDataM,'all')
                    continue
                end
                frameDataS = double(squeeze(allMetricS(:,ifreq,itime)));
                
                
                lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                    [elecIndexM;elecIndexS], 'VariableNames',{'Y','Cond','Sub','Elec'});
                lmeTBL.Cond = nominal(lmeTBL.Cond);
                lmeTBL.Sub = nominal(lmeTBL.Sub);
                lmeTBL.Elec = nominal(lmeTBL.Elec);
                lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
                
                [~,~,lmeStats] = fixedEffects(lmeStruct);
                
                tMap(ifreq,itime) = lmeStats.tStat(2);
                pMap(ifreq,itime) = lmeStats.pValue(2);
                [randBeta,~,~] = randomEffects(lmeStruct);
                Z = designMatrix(lmeStruct,'random');
                Ycorr = lmeTBL.Y-Z*randBeta;
                obsVal1 = Ycorr(lmeTBL.Cond=='1');
                obsVal2 = Ycorr(lmeTBL.Cond=='2');
                y2plot(:,ifreq,itime) = [mean(obsVal1);mean(obsVal2)];
                se2plot(:,ifreq,itime) = [std(obsVal1)./sqrt(numel(obsVal1));std(obsVal2)./sqrt(numel(obsVal2))];
                rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
                rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
                yraw(:,ifreq,itime) = [mean(rawVal1);mean(rawVal2)];
                seraw(:,ifreq,itime) = [std(rawVal1)./sqrt(numel(rawVal1));std(rawVal2)./sqrt(numel(rawVal2))];
                
            end
        end
        
        
        %% section4-2: plot freq coherence
        
        % mid time point of choosen time window
        TimePoint = 0.5;
        
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
        [p_fdr, p_masked] = fdr(pMap(corrFreqInd), 0.05);
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
        
        title([ROIText{iseed} '--' ROIText{isearch} '(' ...
            num2str(TimePoint-0.5*Para.timeWin) '-' ...
            num2str(TimePoint+0.5*Para.timeWin) 's)'])
        xlabel('Frequency (Hz)')
        ylabel('Imaginary Coherence')
        
        
        %% section4-3: plot time coherence
        
        % choosen freq window
        FreqRange = [20 30];
%         FreqRange = [60 90];
        
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
        
        title([ROIText{iseed} '--' ROIText{isearch} '(' ...
            num2str(TimePoint-0.5*Para.timeWin) '-' ...
            num2str(TimePoint+0.5*Para.timeWin) 's)'])
        xlabel('Time relative to camera change (sec)')
        ylabel('Imaginary Coherence')
        
        
        
    end
end
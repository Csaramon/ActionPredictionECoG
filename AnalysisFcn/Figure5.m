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
            
            % calculate pairwise Granger Causality
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
            allGrangerM = [];
            allGrangerS = [];
            for itw = min(timeRange):timeStep:max(timeRange)-timeWin
                
                % choose time window
                cfg = [];
                cfg.latency = [itw itw+timeWin];
                trlDataMtmp = ft_selectdata(cfg,trlDataM);
                
                cfg = [];
                cfg.latency = [itw itw+timeWin];
                trlDataStmp = ft_selectdata(cfg,trlDataS);
                
                % calculate fourier spectrum
                cfg            = [];
                cfg.output     = 'fourier';
                cfg.method     = 'mtmfft';
                %                     cfg.foilim     = [0 250]; % need equidistant frequency bins for granger method
                cfg.tapsmofrq  = 4;
                cfg.keeptrials = 'yes';
                ft_warning off
                freqM    = ft_freqanalysis(cfg, trlDataMtmp);
                
                freqS    = ft_freqanalysis(cfg, trlDataStmp);
                
                % calculate Granger Index in Matched Condition
                grangercfg = [];
                grangercfg.method  = 'granger';
                grangercfg.granger.conditional = 'no';
                grangercfg.granger.sfmethod = 'bivariate';
                grangercfg.channelcmb = {trlDataM.label(seedElec) trlDataM.label(searchElec)};
                %                     grangercfg.channelcmb = sigPair;
                grangerM      = ft_connectivityanalysis(grangercfg, freqM);
                
                for ilabel = 1:numel(grangerM.labelcmb)
                    grangerM.labelcmb{ilabel}(find(grangerM.labelcmb{ilabel}=='['):end)=[];
                end
                fromInd = false(size(grangerM.labelcmb,1),1);
                for is = seedElec'
                    fromInd = fromInd+strcmp(trlDataM.label{is},grangerM.labelcmb(:,1));
                end
                fromInd = logical(fromInd);
                toInd = logical(1-fromInd);
                %                     allGrangerM(:,:,in) = (grangerM.grangerspctrm(fromInd,:)-grangerM.grangerspctrm(toInd,:))./ ...
                %                         (grangerM.grangerspctrm(fromInd,:)+grangerM.grangerspctrm(toInd,:));
                allGrangerM(:,:,in) = grangerM.grangerspctrm(fromInd,:);
                
                % calculate Granger Index in Scambled Condition
                grangerS      = ft_connectivityanalysis(grangercfg, freqS);
                
                for ilabel = 1:numel(grangerS.labelcmb)
                    grangerS.labelcmb{ilabel}(find(grangerS.labelcmb{ilabel}=='['):end)=[];
                end
                fromInd = false(size(grangerS.labelcmb,1),1);
                for is = seedElec'
                    fromInd = fromInd+strcmp(trlDataS.label{is},grangerS.labelcmb(:,1));
                end
                fromInd = logical(fromInd);
                toInd = logical(1-fromInd);
                %                     allGrangerS(:,:,in) = (grangerS.grangerspctrm(fromInd,:)-grangerS.grangerspctrm(toInd,:))./ ...
                %                         (grangerS.grangerspctrm(fromInd,:)+grangerS.grangerspctrm(toInd,:));
                allGrangerS(:,:,in) = grangerS.grangerspctrm(fromInd,:);
                
                timePt(in) = itw+0.5*timeWin;
                in = in+1;
            end
            
            allChanCmb = [allChanCmb;grangerM.labelcmb(fromInd)];
            
            
            % count for pairs of eletrodes
            Npair = Npair + size(grangerM.labelcmb(fromInd),1);
            
            
            % generate index for Subject Electrode and Trial
            metricM = allGrangerM;
            metricS = allGrangerS;
            
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
        

        
        %% plot time frequency Granger Index
        
        tMap(tMap==0)=nan;
        pMap(pMap==0)=nan;
        
        iind = 1:size(pMap,1);
        
        pMap = pMap(iind,:);
        tMap = tMap(iind,:);
        y2plot = y2plot(:,iind,:);
        
        cdat = tMap;
        vdat = Para.freq(iind);
        hdat = Para.timePT;
        
        clim = [min(cdat(:)) max(cdat(:))];
        
        % uncorrected
        %     highlight = pMap< 0.05;
        
        % Bofforoni  correction
        % highlight = pMap< 0.05/numel(pMap);
        
        % FDR  correction
        [p_fdr, p_masked] = fdr(pMap, 0.05);
        highlight = p_masked;
        
        % the significant voxels could be outlined with a black contour
        % plot outline
        hf = figure;
        
        subplot(1,2,1)
        hM = pcolor(hdat,vdat,squeeze(y2plot(1,:,:)));shading interp
        ca = get(gca,'CLim');
        ylabel(['Frequency (Hz)'])
        title('Intact')
        
        [x,y] = meshgrid(hdat, vdat);
        x = interp2(x, 2); % change to 4 for round corners
        y = interp2(y, 2); % change to 4 for round corners
        contourlines = highlight==1;
        
        Ls = bwconncomp(contourlines,4);
        lenths = [];
        for ic = 1:numel(Ls.PixelIdxList)
            [xx,yy] = ind2sub(size(contourlines),Ls.PixelIdxList{ic});
            if Para.timePT(max(yy))-Para.timePT(min(yy)) < 0.095
                contourlines(Ls.PixelIdxList{ic}) = 0;
                
            end
            
        end
        
        contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
        dx = mean(diff(x(1, :))); % remove for round corners
        dy = mean(diff(y(:, 1))); % remove for round corners
        hold on
        contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[238 92 66]/255,'LineWidth',2);
        caxis(ca)
        
        subplot(1,2,2)
        hS = pcolor(hdat,vdat,squeeze(y2plot(2,:,:)));shading interp
        caxis(ca)
        xlabel(['Time (s)'])
        title('Scrambled')
        
        ch = colorbar('peer',gca,'EastOutside');
        set(ch,'position',[0.92 0.12 0.03 0.8],'ticks',[],'fontsize',12, ...
            'ticklabels',{})
        set(ch.Label,'string',['← ' ROIText{iseed} '    From    ' ROIText{isearch} ' →'],...
            'position',[0.9 0],'fontsize',12);
        
        suptitle({'Granger Index'; ...
            [ROIText{iseed} '↔' ROIText{isearch} ' (Nsub=' num2str(numel(unique(lmeTBL.Sub))) ' Nelec=' num2str(numel(unique(lmeTBL.Elec))) ')']});
               
        %% section5-2: plot  non-parametric granger causality across freq
        
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
        hM = shadedErrorBar(x2plot, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255,'linewidth',1,'linestyle','-'});
        hS = shadedErrorBar(x2plot, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255,'linewidth',1,'linestyle','-'});
        
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
        
        title([ROIText{iseed} '--' ROIText{isearch} '(' ...
            num2str(TimePoint-0.5*Para.timeWin) '-' ...
            num2str(TimePoint+0.5*Para.timeWin) 's)'])
        xlabel('Frequency (Hz)')
        ylabel({'GC',['from ' ROIText{iseed} ' to ' ROIText{isearch}]})
        
        
    end
end
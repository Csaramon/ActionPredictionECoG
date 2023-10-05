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


seedIndex = 1;
searchIndex = 2;

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
        allCoordinates = [];
        nelec = 1;
        Para.elecposMNI = [];
        Para.time = [];
        Para.freq = [];
        
        % -------- Subject Level -------- %
        for isub = 1:numel(allsub)
            subname = allsub{isub};
            fprintf(['\n Currently calculating subject: ' subname])
            
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
                cfg.foi          = 2:1:120;
                cfg.tapsmofrq  = 6;
                cfg.keeptrials = 'yes';
                cfg.pad='nextpow2';
                
                
                freqM    = ft_freqanalysis(cfg, trlDataMtmp);
                
                freqS    = ft_freqanalysis(cfg, trlDataStmp);
                
                % calculate COH in Matched Condition
                cfg            = [];
                cfg.method     = 'coh';
                cfg.complex = 'absimag';
                for ielec = 1:numel(seedElec)
                    cfg.channelcmb = {trlDataM.label(seedElec(ielec)) trlDataM.label(searchElec)};
                    COHM             = ft_connectivityanalysis(cfg, freqM);
                    allcohM(ielec,:,in) = mean(COHM.cohspctrm,1);
                    % calculate COH in Scambled Condition
                    cfg.channelcmb = {trlDataS.label(seedElec(ielec)) trlDataS.label(searchElec)};
                    COHS             = ft_connectivityanalysis(cfg, freqS);
                    allcohS(ielec,:,in) = mean(COHS.cohspctrm,1);
                end
                
                timePt(in) = itw+0.5*timeWin;
                in = in+1;
            end
            
            allChanCmb = [allChanCmb;seedElec];
            allCoordinates = [allCoordinates;elecposMNI(seedElec,:)];
            
            allMetricM = cat(1,allMetricM,allcohM);
            allMetricS = cat(1,allMetricS,allcohS);
            
            Para.chanCMB{isub} = allChanCmb;
            Para.timeWin = timeWin;
            Para.timeStep = timeStep;
            Para.timePT = timePt;
            Para.freq = freqM.freq;
        end
        
        %% Plot coherence coordinates correlation
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
        
    end
end
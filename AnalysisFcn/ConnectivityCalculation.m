function varargout = ConnectivityCalculation(calculate)

tic;
if nargin < 1
    calculate = 'Granger';
end

% initialize fieldtrip toolbox
addpath('/data00/Chaoyi/toolbox/tools')
addpath('/data00/Chaoyi/toolbox/fieldtrip-20210418/')
ft_defaults


allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub
basePath = '/data00/Chaoyi/ActionPrediction/';
resultPath = [basePath 'Results/'];



%% -------- ROI Npairtion --------

% Region of Interest include:
ROIIndex = {[1,2],[49,50],[51,52],[53,54],[59,60],[61,62],[63,64],[1],[1]};
ROIAtlas = {'fsAnatomyMacro.nii','fsAnatomyMacro.nii','fsAnatomyMacro.nii','fsAnatomyMacro.nii', ...
    'fsAnatomyMacro.nii','fsAnatomyMacro.nii','fsAnatomyMacro.nii','BA44.nii','PFt.nii'};
ROIText = {'Precentral','SuperiorOccipitalGyrus','MiddleOccipitalGyrus',...
    'InferiorOccipitalGyrus','SuperiorParietalLobe','InferiorParietalLobe','SupraMarginal','BA44','PFt'};

% ROIIndex = {[1008,2008],[1029,2029],[1031,2031],[1022,2022],[1024,2024], ...
%     [1028,2028],[1005,2005],[1011,2011]};
% ROIAtlas = {'aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii', ...
%     'aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii'};
% ROIText = {'InferiorParietalLobe','SuperiorParietalLobe','SupraMarginal','Postcentral','Precentral', ...
%     'SuperiorFrontal','Cuneus','LateralOccipital'};
roiDist = 1; % maximum distance between electrodes and ROI voxels

seedIndex = [1 3 7];
searchIndex = [1 3 7];
allPair = nchoosek(seedIndex,2);
for iseed = seedIndex
    for isearch = searchIndex
        
        %         skip redundant pairs
        if ~ismember([iseed,isearch],allPair,'rows')
            continue
        end
        
        if iseed==isearch
            continue
        end
        
        % load altlas infomation for seed roi
        aparc_nii = load_nifti([basePath 'Atlas' filesep ROIAtlas{iseed}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        if ndims(aparc_vol) == 3
            for iroi = ROIIndex{iseed}
                temp_indx = find(aparc_vol==iroi);
                allroi_index = [allroi_index;temp_indx];
            end
        elseif ndims(aparc_vol) == 4
            %%%------ For probabilistic map ------%%%
            aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{iseed}),4));
            allroi_index = find(aparc_vol>=10);  % minimum probability
        end
        
        % atlas parcellation coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        seed_coordiantes = ras_coordiantes;
        
        % load altlas infomation for search roi
        aparc_nii = load_nifti([basePath 'Atlas' filesep ROIAtlas{isearch}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        if ndims(aparc_vol) == 3
            for iroi = ROIIndex{isearch}
                temp_indx = find(aparc_vol==iroi);
                allroi_index = [allroi_index;temp_indx];
            end
        elseif ndims(aparc_vol) == 4
            %%%------ For probabilistic map ------%%%
            aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{isearch}),4));
            allroi_index = find(aparc_vol>=10);  % minimum probability
        end
        
        % atlas parcellation coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        search_coordiantes = ras_coordiantes;
        
        
        % initialize result variables
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
        
        
        %% -------- Subject Level --------
        
        for isub = 1:numel(allsub)%[1,3,4,6,7,8,9,10]%1:numel(allsub)
            subname = allsub{isub};
            subPath = [basePath subname filesep];
            dataPath = [subPath filesep 'Analysis' filesep];
            fprintf(['\n Currently calculating subject: ' subname])
            
            %% section1: calculate Phase Locking Value %%
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'PLV')
                
                p=0.05; % threshold for IVC
                timeWin = 0.2; % unit in second
                timeStep = 0.05; % unit in second
                freqRange = [3 8];
                
                %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                datafile = dir([resultPath 'TFR' filesep subname filesep 'FreqData.mat']);
                load([datafile.folder filesep datafile.name]);
                
                if exist([dataPath subname 'IVC.mat'],'file')
                    load([dataPath subname 'IVC']);
                end
                
                respElecInd = find(IVC.Intact.theta(:,2)<p | IVC.Intact.alpha(:,2)<p | IVC.Intact.beta(:,2)<p | ...
                    IVC.Intact.lgamma(:,2)<p | IVC.Intact.hgamma(:,2)<p | IVC.Scamble.theta(:,2)<p | ...
                    IVC.Scamble.alpha(:,2)<p | IVC.Scamble.beta(:,2)<p | IVC.Scamble.lgamma(:,2)<p | IVC.Scamble.hgamma(:,2)<p);
                
                
                % choose seed electrodes according to MNI coordinates
                elecposMNI = freqM.elec.elecposMNI;
                tempdev = pdist2(elecposMNI,seed_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                seedElec = unique(it);
                % choose search electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,search_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                searchElec = unique(it);
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec)
                    continue
                end
                
                
                % calculate pairwise PLV
                allPLVM = [];
                allPLVS = [];
                Npair = 1;
                allChanCmb = [];
                timePoint = [];
                timeTFR = freqM.time(2)-freqM.time(1);
                freqInd = find(freqM.freq>=min(freqRange) & freqM.freq<=max(freqRange));
                % calculate connectivity metric
                for iseedElec = seedElec'
                    for isearchElec = searchElec'
                        % ignore electrodes with zero values
                        if ~any(freqM.fourierspctrm(:,iseedElec,:,:)) | ~any(freqM.fourierspctrm(:,isearchElec,:,:))
                            continue
                        end
                        % set all channel combinations
                        %                     cfg.channelcmb = [cfg.channelcmb;[freqM.label(iseedElec),freqM.label(isearchElec)]];
                        
                        % calculate PLV in Matched Condition
                        allChanCmb = [allChanCmb;[iseedElec isearchElec]];
                        
                        timeWinPt = round(1/timeTFR*timeWin);
                        timeStepPt = round(1/timeTFR*timeStep);
                        tmax = floor((size(freqM.time,2)-timeWinPt)/timeStepPt)+1;
                        timePoint = zeros(1,tmax);
                        PLVseriesM = zeros(size(freqM.trialinfo,1),numel(freqInd),tmax);
                        for ifreq = 1:numel(freqInd)
                            for itime = 1:tmax
                                time2cal = round((timeStepPt*(itime-1)+1):(timeStepPt*(itime-1)+timeWinPt));
                                seedFourierM = squeeze(freqM.fourierspctrm(:,iseedElec,freqInd(ifreq),time2cal));
                                searchFourierM = squeeze(freqM.fourierspctrm(:,isearchElec,freqInd(ifreq),time2cal));
                                trlPLV = abs(nanmean(exp(1i*(angle(seedFourierM) -angle(searchFourierM))),2));
                                timePoint(itime) = mean(freqM.time(time2cal));
                                PLVseriesM(:,ifreq,itime) = trlPLV;
                            end
                        end
                        allPLVM(:,Npair,:) = mean(PLVseriesM,2);
                        
                        % calculate PLV in Scrambled Condition
                        PLVseriesS = zeros(size(freqS.trialinfo,1),numel(freqInd),tmax);
                        for ifreq = 1:numel(freqInd)
                            for itime = 1:tmax
                                time2cal = round((timeStepPt*(itime-1)+1):(timeStepPt*(itime-1)+timeWinPt));
                                seedFourierS = squeeze(freqS.fourierspctrm(:,iseedElec,freqInd(ifreq),time2cal));
                                searchFourierS = squeeze(freqS.fourierspctrm(:,isearchElec,freqInd(ifreq),time2cal));
                                trlPLV = abs(nanmean(exp(1i*(angle(seedFourierS) -angle(searchFourierS))),2));
                                PLVseriesS(:,ifreq,itime) = trlPLV;
                            end
                        end
                        allPLVS(:,Npair,:) = mean(PLVseriesS,2);
                        
                        % count for pairs of eletrodes
                        Npair = Npair + 1;
                    end
                end
                
                
                % generate index for Subject Electrode and Trial
                metricM = allPLVM;
                metricS = allPLVS;
                
                elecIndextmp = repmat([nelec:nelec+Npair-2]',1,size(metricM,1))';
                elecIndexM = cat(1,elecIndexM,reshape(elecIndextmp,numel(elecIndextmp),1));
                subIndexM = cat(1,subIndexM,repmat(isub,numel(elecIndextmp),1));
                %                 trlIndexM = cat(1,trlIndexM,repmat(freqM.trialinfo(:,2),size(metricM,2),1));
                trlIndexM = cat(1,trlIndexM,repmat(ones(size(freqM.trialinfo(:,2))),size(metricM,2),1)); % uniform trials
                
                %                 lmeInd = zeros(size(trlDataM.trialinfo(:,2)));
                %                 for im = unique(trlDataM.trialinfo(:,2))'
                %                     a = find(trlDataM.trialinfo(:,2)==im);
                %                     rep = find(a(2:end)-a(1:end-1)>1);
                %                     if isempty(rep)
                %                         lmeInd(a) = a;
                %                     else
                %                         tmpInd = repmat(a(1:rep(1)),numel(rep)+1,1);
                %                         lmeInd(a) = tmpInd(1:numel(a));
                %                     end
                %                 end
                %                 trlIndexM = cat(1,trlIndexM,repmat(lmeInd,size(metricM,2),1));
                
                elecIndextmp = repmat([nelec:nelec+Npair-2]',1,size(metricS,1))';
                elecIndexS = cat(1,elecIndexS,reshape(elecIndextmp,numel(elecIndextmp),1));
                subIndexS = cat(1,subIndexS,repmat(isub,numel(elecIndextmp),1));
                %                 trlIndexS = cat(1,trlIndexS,repmat(freqS.trialinfo(:,2),size(metricS,2),1));
                trlIndexS = cat(1,trlIndexS,repmat(ones(size(freqS.trialinfo(:,2))),size(metricS,2),1)); % uniform trials
                
                %                 lmeInd = zeros(size(trlDataS.trialinfo(:,2)));
                %                 for im = unique(trlDataS.trialinfo(:,2))'
                %                     a = find(trlDataS.trialinfo(:,2)==im);
                %                     rep = find(a(2:end)-a(1:end-1)>1);
                %                     if isempty(rep)
                %                         lmeInd(a) = a;
                %                     else
                %                         tmpInd = repmat(a(1:rep(1)),numel(rep)+1,1);
                %                         lmeInd(a) = tmpInd(1:numel(a));
                %                     end
                %                 end
                %                 trlIndexS = cat(1,trlIndexS,repmat(lmeInd+size(metricM,1),size(metricS,2),1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,reshape(metricM,prod([size(metricM,1),size(metricM,2)]),size(metricM,3)));
                allMetricS = cat(1,allMetricS,reshape(metricS,prod([size(metricS,1),size(metricS,2)]),size(metricS,3)));
                
                Para.chanCMB{isub} = allChanCmb;
                Para.time = timePoint;
                Para.timeTFR = freqM.time;
                Para.freq = freqM.freq;
                nelec = nelec+Npair-1;
                
                
                
            end
            
            %% section1-2: calculate Phase Locking Value across trials%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'PLVtrl')
                
                statP = 0.05; % significance level for statistical test
                clusterP = 0.05; % significance level for cluster
                timeTFR = 0.01; % time step for time frequency results
                
                %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                datafile = dir([resultPath 'TFR' filesep subname filesep 'FreqData.mat']);
                load([datafile.folder filesep datafile.name]);
                
                % choose seed electrodes according to MNI coordinates
                elecposMNI = freqM.elec.elecposMNI;
                tempdev = pdist2(elecposMNI,seed_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                seedElec = unique(it);
                % choose search electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,search_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                searchElec = unique(it);
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec)
                    continue
                end
                
                
                % calculate pairwise PLV
                allPLVM = [];
                allPLVS = [];
                Npair = 1;
                allChanCmb = [];
                % calculate connectivity metric
                for iseedElec = seedElec'
                    for isearchElec = searchElec'
                        % ignore electrodes with zero values
                        if ~any(freqM.fourierspctrm(:,iseedElec,:,:)) | ~any(freqM.fourierspctrm(:,isearchElec,:,:))
                            continue
                        end
                        
                        % calculate PLV in Matched Condition
                        allChanCmb = [allChanCmb;[iseedElec isearchElec]];
                        PLVseriesM = zeros(numel(freqM.freq),numel(freqM.time));
                        for ifreq = 1:numel(freqM.freq)
                            seedFourierM = squeeze(freqM.fourierspctrm(:,iseedElec,ifreq,:))';
                            searchFourierM = squeeze(freqM.fourierspctrm(:,isearchElec,ifreq,:))';
                            PLVM = abs(nanmean(exp(1i*(angle(seedFourierM) -angle(searchFourierM))),2));
                            PLVseriesM(ifreq,:) = PLVM;
                        end
                        allPLVM(Npair,:,:) = PLVseriesM;
                        
                        % calculate PLV in Scrambled Condition
                        PLVseriesS = zeros(numel(freqS.freq),numel(freqS.time));
                        for ifreq = 1:numel(freqS.freq)
                            seedFourierS = squeeze(freqS.fourierspctrm(:,iseedElec,ifreq,:))';
                            searchFourierS = squeeze(freqS.fourierspctrm(:,isearchElec,ifreq,:))';
                            PLVS = abs(nanmean(exp(1i*(angle(seedFourierS) -angle(searchFourierS))),2));
                            PLVseriesS(ifreq,:) = PLVS;
                        end
                        allPLVS(Npair,:,:) = PLVseriesS;
                        
                        % count for pairs of eletrodes
                        Npair = Npair + 1;
                    end
                end
                
                
                % generate index for Subject Electrode and Trial
                metricM = allPLVM;
                metricS = allPLVS;
                
                elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                
                elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,metricM);
                allMetricS = cat(1,allMetricS,metricS);
                
                Para.chanCMB{isub} = allChanCmb;
                Para.time = freqM.time;
                Para.freq = freqM.freq;
                nelec = nelec+Npair-1;
                
                
                
            end
            
            
            %% section2: calculate Coherence %%
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'COH')
                
                p=0.05; % threshold for IVC
                timeWin = [0 1]; % unit in second
                
                %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                datafile = dir([dataPath subname 'LAR_trlData.mat']);
                load([datafile.folder filesep datafile.name]);
                
                if exist([dataPath subname 'IVC.mat'],'file')
                    load([dataPath subname 'IVC']);
                end
                
                respElecInd = find(IVC.Intact.theta(:,2)<p | IVC.Intact.alpha(:,2)<p | IVC.Intact.beta(:,2)<p | ...
                    IVC.Intact.lgamma(:,2)<p | IVC.Intact.hgamma(:,2)<p | IVC.Scamble.theta(:,2)<p | ...
                    IVC.Scamble.alpha(:,2)<p | IVC.Scamble.beta(:,2)<p | IVC.Scamble.lgamma(:,2)<p | IVC.Scamble.hgamma(:,2)<p);
                
                % choose seed electrodes according to MNI coordinates
                elecposMNI = trlData.elec.elecposMNI;
                tempdev = pdist2(elecposMNI,seed_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                seedElec = unique(it);
                %                 seedElec = intersect(seedElec,respElecInd);
                
                % choose search electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,search_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                searchElec = unique(it);
                %                 searchElec = intersect(searchElec,respElecInd);
                
                cfg = [];
                cfg.demean = 'yes';
                cfg.detrend = 'yes';
                
                trlData = ft_preprocessing(cfg,trlData);
                
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
                cfg.latency = timeWin;
                cfg.trials = find(trlData.trialinfo(:,1)==0);
                trlDataM = ft_selectdata(cfg,trlData);
                
                cfg = [];
                cfg.latency = timeWin;
                cfg.trials = find(trlData.trialinfo(:,1)==1);
                trlDataS = ft_selectdata(cfg,trlData);
                
                % calculate fourier spectrum
                cfg            = [];
                cfg.output     = 'fourier';
                cfg.method     = 'mtmfft';
                cfg.foilim     = [2 120];
                % cfg.foi          = logspace(log10(2),log10(128),32);
                cfg.tapsmofrq  = 5;
                cfg.keeptrials = 'yes';
                freqM    = ft_freqanalysis(cfg, trlDataM);
                
                cfg            = [];
                cfg.output     = 'fourier';
                cfg.method     = 'mtmfft';
                cfg.foilim     = [2 120];
                % cfg.foi          = logspace(log10(2),log10(128),32);
                cfg.tapsmofrq  = 5;
                cfg.keeptrials = 'yes';
                freqS    = ft_freqanalysis(cfg, trlDataS);
                
                
                
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
                % calculate COH in Matched Condition
                
                cfg            = [];
                cfg.method     = 'coh';
                %                 cfg.complex = 'absimag';
                cfg.channelcmb = {trlData.label(seedElec) trlData.label(searchElec)};
                COHM             = ft_connectivityanalysis(cfg, freqM);
                
                cfg            = [];
                cfg.method     = 'coh';
                %                 cfg.complex = 'absimag';
                cfg.channelcmb = {trlData.label(seedElec) trlData.label(searchElec)};
                COHS             = ft_connectivityanalysis(cfg, freqS);
                
                allChanCmb = [allChanCmb;COHM.labelcmb];
                
                
                % count for pairs of eletrodes
                Npair = Npair + size(COHM.labelcmb,1);
                
                
                
                % generate index for Subject Electrode and Trial
                metricM = COHM.cohspctrm;
                metricS = COHS.cohspctrm;
                
                elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                
                elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,metricM);
                allMetricS = cat(1,allMetricS,metricS);
                
                
                Para.chanCMB{isub} = allChanCmb;
                Para.timeWin = timeWin;
                Para.freq = freqM.freq;
                nelec = nelec+Npair-1;
                
                
                
            end
            
            
            %% section2-2: calculate Coherence across trials%%
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'COHtrl')
                
                p=0.05; % threshold for IVC
                
                %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                datafile = dir([resultPath 'TFR' filesep subname filesep 'FreqData.mat']);
                load([datafile.folder filesep datafile.name]);
                
                if exist([dataPath subname 'IVC.mat'],'file')
                    load([dataPath subname 'IVC']);
                end
                
                respElecInd = find(IVC.Intact.theta(:,2)<p | IVC.Intact.alpha(:,2)<p | IVC.Intact.beta(:,2)<p | ...
                    IVC.Intact.lgamma(:,2)<p | IVC.Intact.hgamma(:,2)<p | IVC.Scamble.theta(:,2)<p | ...
                    IVC.Scamble.alpha(:,2)<p | IVC.Scamble.beta(:,2)<p | IVC.Scamble.lgamma(:,2)<p | IVC.Scamble.hgamma(:,2)<p);
                
                % choose seed electrodes according to MNI coordinates
                elecposMNI = freqM.elec.elecposMNI;
                tempdev = pdist2(elecposMNI,seed_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                seedElec = unique(it);
                seedElec = intersect(seedElec,respElecInd);
                
                % choose search electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,search_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                searchElec = unique(it);
                searchElec = intersect(searchElec,respElecInd);
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec)
                    continue
                end
                
                
                % calculate pairwise COH
                allCOHM = [];
                allCOHS = [];
                Npair = 1;
                allChanCmb = [];
                timePoint = [];
                % calculate connectivity metric
                for iseedElec = seedElec'
                    for isearchElec = searchElec'
                        % ignore electrodes with zero values
                        if ~any(freqM.fourierspctrm(:,iseedElec,:,:)) | ~any(freqM.fourierspctrm(:,isearchElec,:,:))
                            continue
                        end
                        
                        % calculate COH in Matched Condition
                        allChanCmb = [allChanCmb;[iseedElec isearchElec]];
                        COHseriesM = zeros(numel(freqM.freq),numel(freqM.time));
                        for ifreq = 1:numel(freqM.freq)
                            seedFourierM = squeeze(freqM.fourierspctrm(:,iseedElec,ifreq,:))';
                            searchFourierM = squeeze(freqM.fourierspctrm(:,isearchElec,ifreq,:))';
                            spec1 = nanmean(seedFourierM.*conj(seedFourierM),2);
                            spec2 = nanmean(searchFourierM.*conj(searchFourierM),2);
                            specX = abs(nanmean(seedFourierM.*conj(searchFourierM),2)).^2;
                            COHseriesM(ifreq,:) = specX./(spec1.*spec2);
                        end
                        allCOHM(Npair,:,:) = COHseriesM;
                        
                        % calculate COH in Scrambled Condition
                        COHseriesS = zeros(numel(freqS.freq),numel(freqS.time));
                        for ifreq = 1:numel(freqS.freq)
                            seedFourierS = squeeze(freqS.fourierspctrm(:,iseedElec,ifreq,:))';
                            searchFourierS = squeeze(freqS.fourierspctrm(:,isearchElec,ifreq,:))';
                            spec1 = nanmean(seedFourierS.*conj(seedFourierS),2);
                            spec2 = nanmean(searchFourierS.*conj(searchFourierS),2);
                            specX = abs(nanmean(seedFourierS.*conj(searchFourierS),2)).^2;
                            COHseriesS(ifreq,:) = specX./(spec1.*spec2);
                        end
                        allCOHS(Npair,:,:) = COHseriesS;
                        
                        % count for pairs of eletrodes
                        Npair = Npair + 1;
                    end
                end
                
                
                % generate index for Subject Electrode and Trial
                metricM = allCOHM;
                metricS = allCOHS;
                
                elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                
                elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,metricM);
                allMetricS = cat(1,allMetricS,metricS);
                
                Para.chanCMB{isub} = allChanCmb;
                Para.time = freqM.time;
                Para.freq = freqM.freq;
                nelec = nelec+Npair-1;
                
                
                
            end
            
            
            
            %% section3: calculate Phase Slope Index %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'PSI')
                
                timeTFR = 0.002; % time step for time frequency results
                timeWin = 0.2; % unit in second
                timeStep = 0.05; % unit in second
                freqRange = [60 120];
                
                %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if exist([dataPath subname 'LAR_trlData.mat'],'file')
                    a = load([dataPath subname 'LAR_trlData']);
                end
                c = fieldnames(a);
                trlData = a.(c{1});
                
                % choose seed electrodes according to MNI coordinates
                elecposMNI = trlData.elec.elecposMNI;
                tempdev = pdist2(elecposMNI,seed_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                seedElec = unique(it);
                % choose search electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,search_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                searchElec = unique(it);
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec)
                    continue
                end
                
                cfg = [];
                cfg.demean = 'yes';
                cfg.detrend = 'yes';
                
                trlData = ft_preprocessing(cfg,trlData);
                
                badChanInd = trlData.trial{1,1}(seedElec,1)==0;
                seedElec(badChanInd) = [];
                
                badChanInd = trlData.trial{1,1}(searchElec,1)==0;
                searchElec(badChanInd) = [];
                
                % remove duplicated elec index
                duplicateInd = intersect(searchElec,seedElec);
                if ~isempty(duplicateInd)
                    if numel(searchElec) > numel(seedElec)
                        searchElec = setdiff(searchElec,duplicateInd);
                    else
                        seedElec = setdiff(seedElec,duplicateInd);
                    end
                end
                % seperate conditions
                cfg = [];
                %                 cfg.latency = [-0.5 1];
                cfg.trials = find(trlData.trialinfo(:,1)==0);
                cfg.channel = [seedElec;searchElec];
                trlDataM = ft_selectdata(cfg,trlData);
                
                cfg = [];
                %                 cfg.latency = [-0.5 1];
                cfg.trials = find(trlData.trialinfo(:,1)==1);
                cfg.channel = [seedElec;searchElec];
                trlDataS = ft_selectdata(cfg,trlData);
                clear trlData
                
                
                %%%---------- hilbert ----------%%%
                cfg              = [];
                cfg.output       = 'fourier';
                cfg.method       = 'hilbert';
                cfg.foi          = freqRange(1):5:freqRange(2);
                cfg.width        =  5;
                cfg.toi          = min(trlDataM.time{1}):timeTFR:max(trlDataM.time{1}); %'all';
                cfg.keeptrials   = 'yes';
                cfg.filttype = 'firws';
                cfg.filtorder = nan;
                cfg.filtdir = 'onepass-zerophase';
                cfg.precision = 'single';
                
                ft_warning off
                freqM = ft_freqanalysis(cfg,trlDataM);
                freqS = ft_freqanalysis(cfg,trlDataS);
                
                % calculate pairwise PSI
                allPSIM = [];
                allPSIS = [];
                Npair = 1;
                allChanCmb = [];
                timePoint = [];
                % calculate connectivity metric
                for iseedElec = 1:numel(seedElec)
                    for isearchElec = numel(seedElec)+1:numel(seedElec)+numel(searchElec)
                        
                        % calculate PSI in Matched Condition
                        allElecInd = [seedElec;searchElec];
                        allChanCmb = [allChanCmb;[allElecInd(iseedElec) allElecInd(isearchElec)]];
                        timeWinPt = round(1/timeTFR*timeWin);
                        timeStepPt = round(1/timeTFR*timeStep);
                        tmax = floor((size(freqM.time,2)-timeWinPt)/timeStepPt)+1;
                        timePoint = zeros(1,tmax);
                        PSIseriesM = zeros(1,tmax);
                        
                        for itime = 1:tmax
                            time2cal = round((timeStepPt*(itime-1)+1):(timeStepPt*(itime-1)+timeWinPt));
                            seedFourierM = squeeze(nanmean(freqM.fourierspctrm(:,iseedElec,:,time2cal),4));
                            searchFourierM = squeeze(nanmean(freqM.fourierspctrm(:,isearchElec,:,time2cal),4));
                            cs  = zeros(2,2,size(seedFourierM,2));
                            for itrl = 1:size(seedFourierM,1)
                                for ifreq = 1:size(seedFourierM,2)
                                    
                                    tempfftdat = [seedFourierM(itrl,ifreq),searchFourierM(itrl,ifreq)];
                                    cs(:,:,ifreq) = cs(:,:,ifreq) + tempfftdat'*tempfftdat;
                                    
                                end
                            end
                            cs = cs./itrl;
                            
                            pp = zeros(2,2,size(seedFourierM,2));
                            
                            for fi=1:size(seedFourierM,2)
                                pp(:,:,fi) = cs(:,:,fi)./sqrt(diag(cs(:,:,fi))*diag(cs(:,:,fi))');
                            end
                            
                            psi_observe = sum(imag(conj(pp(:,:,1:end-1)).*pp(:,:,2:end)),3);
                            PSIseriesM(1,itime) = psi_observe(1,2);
                            timePoint(itime) = mean(freqM.time(time2cal));
                        end
                        allPSIM(Npair,:) = PSIseriesM;
                        
                        % calculate COH in Scrambled Condition
                        PSIseriesS = zeros(1,tmax);
                        for itime = 1:tmax
                            time2cal = round((timeStepPt*(itime-1)+1):(timeStepPt*(itime-1)+timeWinPt));
                            seedFourierS = squeeze(nanmean(freqS.fourierspctrm(:,iseedElec,:,time2cal),4));
                            searchFourierS = squeeze(nanmean(freqS.fourierspctrm(:,isearchElec,:,time2cal),4));
                            cs  = zeros(2,2,size(seedFourierS,2));
                            for itrl = 1:size(seedFourierS,1)
                                for ifreq = 1:size(seedFourierS,2)
                                    
                                    tempfftdat = [seedFourierS(itrl,ifreq),searchFourierS(itrl,ifreq)];
                                    cs(:,:,ifreq) = cs(:,:,ifreq) + tempfftdat'*tempfftdat;
                                    
                                end
                            end
                            cs = cs./itrl;
                            
                            pp = zeros(2,2,size(seedFourierS,2));
                            
                            for fi=1:size(seedFourierS,2)
                                pp(:,:,fi) = cs(:,:,fi)./sqrt(diag(cs(:,:,fi))*diag(cs(:,:,fi))');
                            end
                            
                            psi_observe = sum(imag(conj(pp(:,:,1:end-1)).*pp(:,:,2:end)),3);
                            PSIseriesS(1,itime) = psi_observe(1,2);
                        end
                        allPSIS(Npair,:) = PSIseriesS;
                        
                        % count for pairs of eletrodes
                        Npair = Npair + 1;
                    end
                end
                
                
                % generate index for Subject Electrode and Trial
                metricM = allPSIM;
                metricS = allPSIS;
                
                elecIndexM =  cat(1,elecIndexM,isub*1000+allChanCmb);
                %                 elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                elecIndexS = cat(1,elecIndexS,isub*1000+allChanCmb);
                %                 elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,metricM);
                allMetricS = cat(1,allMetricS,metricS);
                
                Para.chanCMB{isub} = allChanCmb;
                Para.time = timePoint;
                Para.timeTFR = timeTFR;
                Para.freq = freqRange;
                nelec = nelec+Npair-1;
                
                
                
            end
            
            %% section4: calculate multivariate granger causality (MVGC) %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(calculate,'mvgc')
                
                timeWin = [0 1]; % unit in second
                
                %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if exist([dataPath subname 'LAR_trlData.mat'],'file')
                    a = load([dataPath subname 'LAR_trlData']);
                end
                c = fieldnames(a);
                trlData = a.(c{1});
                
                % choose seed electrodes according to MNI coordinates
                elecposMNI = trlData.elec.elecposMNI;
                tempdev = pdist2(elecposMNI,seed_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                seedElec = unique(it);
                % choose search electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,search_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                searchElec = unique(it);
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec)
                    continue
                end
                
                cfg = [];
                cfg.demean = 'yes';
                cfg.detrend = 'yes';
                
                trlData = ft_preprocessing(cfg,trlData);
                
                badChanInd = trlData.trial{1,1}(seedElec,1)==0;
                seedElec(badChanInd) = [];
                
                badChanInd = trlData.trial{1,1}(searchElec,1)==0;
                searchElec(badChanInd) = [];
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec) | isempty(setdiff(seedElec,searchElec))
                    continue
                end
                
                % remove superimposed electrodes
                duplicateInd = intersect(searchElec,seedElec);
                if ~isempty(duplicateInd)
                    if numel(searchElec) > numel(seedElec)
                        searchElec = setdiff(searchElec,duplicateInd);
                    else
                        seedElec = setdiff(seedElec,duplicateInd);
                    end
                end
                % seperate conditions
                cfg = [];
                cfg.trials = find(trlData.trialinfo(:,1)==0);
                cfg.channel = [seedElec;searchElec];
                trlDataM = ft_selectdata(cfg,trlData);
                
                cfg = [];
                cfg.trials = find(trlData.trialinfo(:,1)==1);
                cfg.channel = [seedElec;searchElec];
                trlDataS = ft_selectdata(cfg,trlData);
                clear trlData
                
                
                % calculate mvgc
                allmvgcM = [];
                allmvgcS = [];
                Npair = 1;
                allChanCmb = [];
                
                time2cal  = trlDataM.time{1}>=min(timeWin) & trlDataM.time{1}<=max(timeWin);
                
                X = [];
                for itrl = 1:numel(trlDataM.trial)
                    X(:,:,itrl) = trlDataM.trial{itrl}(:,time2cal);
                end
                
                ntrials   = size(X,3);     % number of trials
                nobs      = size(X,2);   % number of observations per trial
                mvgc_calculate % calculate granger causality
                FM = f;
                fres = size(f,3)-1;
                freqPoints = sfreqs(fres,fs)';
                morderM = morder+1;
                
                X = [];
                for itrl = 1:numel(trlDataS.trial)
                    X(:,:,itrl) = trlDataS.trial{itrl}(:,time2cal);
                end
                
                ntrials   = size(X,3);     % number of trials
                nobs      = size(X,2);   % number of observations per trial
                mvgc_calculate % calculate granger causality
                FS = f;
                morderS = morder+1;
                % calculate connectivity metric
                for iseedElec = 1:numel(seedElec)
                    for isearchElec = numel(seedElec)+1:numel(seedElec)+numel(searchElec)
                        % ignore electrodes with zero values
                        if isnan(FM(isearchElec,iseedElec,:)) | isnan(FS(isearchElec,iseedElec,:)) ...
                                | ~any(FM(isearchElec,iseedElec,:)) | ~any(FS(isearchElec,iseedElec,:))
                            continue
                        end
                        allElecInd = [seedElec;searchElec];
                        allChanCmb = [allChanCmb;[allElecInd(iseedElec) allElecInd(isearchElec)]];
                        allmvgcM(Npair,:) = FM(isearchElec,iseedElec,:);
                        
                        allmvgcS(Npair,:) = FS(isearchElec,iseedElec,:);
                        
                        % count for pairs of eletrodes
                        Npair = Npair + 1;
                    end
                end
                
                
                % generate index for Subject Electrode and Trial
                metricM = allmvgcM;
                metricS = allmvgcS;
                
                elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                
                elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,metricM);
                allMetricS = cat(1,allMetricS,metricS);
                
                Para.chanCMB{isub} = allChanCmb;
                Para.morder{isub} = [morderM,morderS];
                Para.freq = freqPoints;
                nelec = nelec+Npair-1;
                
                
            end
            
            
            %% section4-2: calculate Granger causality %%
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'Granger')
                
                p=0.05; % threshold for IVC
                timeWin = [0 0.5]; % unit in second
                
                %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                datafile = dir([dataPath subname 'LAR_trlData.mat']);
                load([datafile.folder filesep datafile.name]);
                
                if exist([dataPath subname 'IVC.mat'],'file')
                    load([dataPath subname 'IVC']);
                end
                
                respElecInd = find(IVC.Intact.theta(:,2)<p | IVC.Intact.alpha(:,2)<p | IVC.Intact.beta(:,2)<p | ...
                    IVC.Intact.lgamma(:,2)<p | IVC.Intact.hgamma(:,2)<p | IVC.Scamble.theta(:,2)<p | ...
                    IVC.Scamble.alpha(:,2)<p | IVC.Scamble.beta(:,2)<p | IVC.Scamble.lgamma(:,2)<p | IVC.Scamble.hgamma(:,2)<p);
                
                % choose seed electrodes according to MNI coordinates
                elecposMNI = trlData.elec.elecposMNI;
                tempdev = pdist2(elecposMNI,seed_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                seedElec = unique(it);
                %                 seedElec = intersect(seedElec,respElecInd);
                
                % choose search electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,search_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                searchElec = unique(it);
                %                 searchElec = intersect(searchElec,respElecInd);
                
                cfg = [];
                cfg.demean = 'yes';
                cfg.detrend = 'yes';
                
                trlData = ft_preprocessing(cfg,trlData);
                
                % skip bad channels
                badChanInd = trlData.trial{1,1}(seedElec,1)==0;
                seedElec(badChanInd) = [];
                
                badChanInd = trlData.trial{1,1}(searchElec,1)==0;
                searchElec(badChanInd) = [];
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec) | isempty(setdiff(seedElec,searchElec))
                    continue
                end
                
                % remove superimposed electrodes
                duplicateInd = intersect(searchElec,seedElec);
                if ~isempty(duplicateInd)
                    if numel(searchElec) > numel(seedElec)
                        searchElec = setdiff(searchElec,duplicateInd);
                    else
                        seedElec = setdiff(seedElec,duplicateInd);
                    end
                end
                % seperate conditions
                cfg = [];
                cfg.latency = timeWin;
                cfg.trials = find(trlData.trialinfo(:,1)==0);
                trlDataM = ft_selectdata(cfg,trlData);
                
                cfg = [];
                cfg.latency = timeWin;
                cfg.trials = find(trlData.trialinfo(:,1)==1);
                trlDataS = ft_selectdata(cfg,trlData);
                
                
                % calculate pairwise Granger Index
                Npair = 1;
                allChanCmb = [];
                allGrangerM = [];
                allGrangerS = [];
                for iseedElec = seedElec'
                    for isearchElec = searchElec'
                        allChanCmb = [allChanCmb;[iseedElec,isearchElec]];
                        % Matched Condition
                        
                        % calculate multivariate autoregressive model
                        cfg         = [];
                        cfg.order   = 10;
                        cfg.toolbox = 'bsmart';
                        cfg.channel = [iseedElec,isearchElec];
                        mdata       = ft_mvaranalysis(cfg, trlDataM);
                        % calculate spectral transfer function
                        cfg        = [];
                        cfg.method = 'mvar';
                        mfreq      = ft_freqanalysis(cfg, mdata);
                        % calculate granger index
                        cfg            = [];
                        cfg.method     = 'granger';
                        grangerM             = ft_connectivityanalysis(cfg, mfreq);
                        if strcmp(trlData.label(iseedElec),grangerM.label(1))
                            allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(1,2,:)-grangerM.grangerspctrm(2,1,:))./ ...
                                (grangerM.grangerspctrm(1,2,:)+grangerM.grangerspctrm(2,1,:)))';
                            
                        elseif strcmp(trlData.label(iseedElec),grangerM.label(2))
                            allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(2,1,:)-grangerM.grangerspctrm(1,2,:))./ ...
                                (grangerM.grangerspctrm(1,2,:)+grangerM.grangerspctrm(2,1,:)))';
                        end
                        
                        %                         cfg = [];
                        %                         cfg.method = 'mtmfft';
                        %                         cfg.output = 'fourier';
                        %                         cfg.channel = [iseedElec,isearchElec];
                        %                         cfg.tapsmofrq = 2;
                        %                         freqdataM = ft_freqanalysis(cfg, trlDataM);
                        %
                        %                         grangercfg = [];
                        %                         grangercfg.method  = 'granger';
                        %                         grangercfg.granger.conditional = 'no';
                        %                         grangercfg.granger.sfmethod = 'bivariate';
                        %
                        %                         grangerM      = ft_connectivityanalysis(grangercfg, freqdataM);
                        %
                        %                         if strncmp(trlData.label{iseedElec},grangerM.labelcmb(1),length(trlData.label{iseedElec}))
                        %                             allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(1,:)-grangerM.grangerspctrm(2,:))./ ...
                        %                                 (grangerM.grangerspctrm(1,:)+grangerM.grangerspctrm(2,:)))';
                        %                         elseif strncmp(trlData.label{iseedElec},grangerM.labelcmb(2),length(trlData.label{iseedElec}))
                        %                             allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(2,:)-grangerM.grangerspctrm(1,:))./ ...
                        %                                 (grangerM.grangerspctrm(1,:)+grangerM.grangerspctrm(2,:)))';
                        %                         end
                        
                        
                        % Scrambled Condition
                        
                        % calculate multivariate autoregressive model
                        cfg         = [];
                        cfg.order   = 10;
                        cfg.toolbox = 'bsmart';
                        cfg.channel = [iseedElec,isearchElec];
                        mdata       = ft_mvaranalysis(cfg, trlDataS);
                        % calculate spectral transfer function
                        cfg        = [];
                        cfg.method = 'mvar';
                        mfreq      = ft_freqanalysis(cfg, mdata);
                        % calculate granger index
                        cfg            = [];
                        cfg.method     = 'granger';
                        grangerS             = ft_connectivityanalysis(cfg, mfreq);
                        if strcmp(trlData.label(iseedElec),grangerS.label(1))
                            allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(1,2,:)-grangerS.grangerspctrm(2,1,:))./ ...
                                (grangerS.grangerspctrm(1,2,:)+grangerS.grangerspctrm(2,1,:)))';
                            
                        elseif strcmp(trlData.label(iseedElec),grangerS.label(2))
                            allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(2,1,:)-grangerS.grangerspctrm(1,2,:))./ ...
                                (grangerS.grangerspctrm(1,2,:)+grangerS.grangerspctrm(2,1,:)))';
                        end
                        %                         cfg = [];
                        %                         cfg.method = 'mtmfft';
                        %                         cfg.output = 'fourier';
                        %                         cfg.channel = [iseedElec,isearchElec];
                        %                         cfg.tapsmofrq = 2;
                        %                         freqdataS = ft_freqanalysis(cfg, trlDataS);
                        %
                        %                         grangercfg = [];
                        %                         grangercfg.method  = 'granger';
                        %                         grangercfg.granger.conditional = 'no';
                        %                         grangercfg.granger.sfmethod = 'bivariate';
                        %
                        %                         grangerS      = ft_connectivityanalysis(grangercfg, freqdataS);
                        %
                        %                         if strncmp(trlData.label{iseedElec},grangerS.labelcmb(1),length(trlData.label{iseedElec}))
                        %                             allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(1,:)-grangerS.grangerspctrm(2,:))./ ...
                        %                                 (grangerS.grangerspctrm(1,:)+grangerS.grangerspctrm(2,:)))';
                        %                         elseif strncmp(trlData.label{iseedElec},grangerS.labelcmb(2),length(trlData.label{iseedElec}))
                        %                             allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(2,:)-grangerS.grangerspctrm(1,:))./ ...
                        %                                 (grangerS.grangerspctrm(1,:)+grangerS.grangerspctrm(2,:)))';
                        %                         end
                        
                        
                        % count for pairs of eletrodes
                        Npair = Npair + 1;
                        
                        
                    end
                end
                
                % generate index for Subject Electrode and Trial
                
                elecIndexM =  cat(1,elecIndexM,isub*1000+allChanCmb);
                %                 elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                elecIndexS = cat(1,elecIndexS,isub*1000+allChanCmb);
                %                 elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,allGrangerM);
                allMetricS = cat(1,allMetricS,allGrangerS);
                
                
                Para.chanCMB{isub} = allChanCmb;
                Para.timeWin = timeWin;
                Para.freq = grangerM.freq;
                nelec = nelec+Npair-1;
                
                
                
            end
            
            
            
        end
        
        
        
        
        %% -------- Group Level --------
        
        %% section1: calculate Phase Locking Value (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PLV')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(1,size(metricM,3));
            pMap = zeros(1,size(metricM,3));
            y2plot = zeros(2,size(metricM,3));
            se2plot = zeros(2,size(metricM,3));
            
            strlen = 0;
            % TFR point wise LME
            for itime = 1:size(metricM,3)
                
                
                s = ['Calculating tf point: ' num2str(itime) '/' num2str(size(metricM,3)) 'times'];
                strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
                strlen = strlentmp - strlen;
                
                frameDataM = double(squeeze(allMetricM(:,itime)));
                % skip nan point
                if ~any(frameDataM,'all')
                    continue
                end
                frameDataS = double(squeeze(allMetricS(:,itime)));
                
                
                lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                    [elecIndexM;elecIndexS],[trlIndexM;trlIndexS], 'VariableNames',{'Y','Cond','Sub','Elec','Trl'});
                lmeTBL.Cond = nominal(lmeTBL.Cond);
                lmeTBL.Sub = nominal(lmeTBL.Sub);
                lmeTBL.Elec = nominal(lmeTBL.Elec);
                lmeTBL.Trl = nominal(lmeTBL.Trl);
                lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)+(1|Trl)','fitmethod','reml','DummyVarCoding','effects');
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
                
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz' ...
                filesep ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','Para');
            
            
        end
        
        %% section1-2: calculate Phase Locking Value across trials(LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PLVtrl')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(size(metricM,2),size(metricM,3));
            pMap = zeros(size(metricM,2),size(metricM,3));
            
            strlen = 0;
            % TFR point wise LME
            for ifreq = 1:size(metricM,2)
                for itime = 1:size(metricM,3)
                    
                    s = ['Calculating tf point: ' num2str(ifreq) '/' num2str(size(metricM,2)) 'freqs, ' ...
                        num2str(itime) '/' num2str(size(metricM,3)) 'times'];
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
                    
                end
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','Para');
            
            
            
        end
        
        
        %% section2: calculate Coherence (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'COH')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(1,size(metricM,2));
            pMap = zeros(1,size(metricM,2));
            y2plot = zeros(2,size(metricM,2));
            se2plot = zeros(2,size(metricM,2));
            
            strlen = 0;
            % TFR point wise LME
            for ifreq = 1:size(metricM,2)
                
                
                s = ['Calculating freq point: ' num2str(ifreq) '/' num2str(size(metricM,2)) 'freqs'];
                strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
                strlen = strlentmp - strlen;
                
                frameDataM = double(squeeze(allMetricM(:,ifreq)));
                % skip nan point
                if ~any(frameDataM,'all')
                    continue
                end
                frameDataS = double(squeeze(allMetricS(:,ifreq)));
                
                
                lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                    [elecIndexM;elecIndexS], 'VariableNames',{'Y','Cond','Sub','Elec'});
                lmeTBL.Cond = nominal(lmeTBL.Cond);
                lmeTBL.Sub = nominal(lmeTBL.Sub);
                lmeTBL.Elec = nominal(lmeTBL.Elec);
                lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
                [~,~,lmeStats] = fixedEffects(lmeStruct);
                tMap(1,ifreq) = lmeStats.tStat(2);
                pMap(1,ifreq) = lmeStats.pValue(2);
                [randBeta,~,~] = randomEffects(lmeStruct);
                Z = designMatrix(lmeStruct,'random');
                Ycorr = lmeTBL.Y-Z*randBeta;
                obsVal1 = Ycorr(lmeTBL.Cond=='1');
                obsVal2 = Ycorr(lmeTBL.Cond=='2');
                y2plot(:,ifreq) = [mean(obsVal1);mean(obsVal2)];
                se2plot(:,ifreq) = [std(obsVal1)./sqrt(numel(obsVal1));std(obsVal2)./sqrt(numel(obsVal2))];
                
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','Para');
            
            
        end
        
        
        %% section2-2: calculate Coherence across trials(LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'COHtrl')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(size(metricM,2),size(metricM,3));
            pMap = zeros(size(metricM,2),size(metricM,3));
            y2plot = zeros(2,size(metricM,2),size(metricM,3));
            se2plot = zeros(2,size(metricM,2),size(metricM,3));
            
            strlen = 0;
            % TFR point wise LME
            for ifreq = 1:size(metricM,2)
                for itime = 1:size(metricM,3)
                    
                    s = ['Calculating tf point: ' num2str(ifreq) '/' num2str(size(metricM,2)) 'freqs, ' ...
                        num2str(itime) '/' num2str(size(metricM,3)) 'times'];
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
                    
                end
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','Para');
            
            
        end
        
        
        %% section3: calculate Phase Slope Index (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PSI')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(1,size(metricM,3));
            pMap = zeros(1,size(metricM,3));
            y2plot = zeros(2,size(metricM,3));
            y2raw = zeros(2,size(metricM,3));
            se2plot = zeros(2,size(metricM,3));
            
            strlen = 0;
            % PSI point wise LME
            for itime = 1:size(metricM,2)
                
                
                s = ['Calculating tf point: ' num2str(itime) '/' num2str(size(metricM,2)) 'times'];
                strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
                strlen = strlentmp - strlen;
                
                frameDataM = double(squeeze(allMetricM(:,itime)));
                % skip nan point
                if ~any(frameDataM,'all')
                    continue
                end
                frameDataS = double(squeeze(allMetricS(:,itime)));
                
                
                %                 lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                %                     [elecIndexM;elecIndexS], 'VariableNames',{'Y','Cond','Sub','Elec'});
                %                 lmeTBL.Cond = nominal(lmeTBL.Cond);
                %                 lmeTBL.Sub = nominal(lmeTBL.Sub);
                %                 lmeTBL.Elec = nominal(lmeTBL.Elec);
                %                 lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
                
                lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                    [elecIndexM(:,1);elecIndexS(:,1)], [elecIndexM(:,2);elecIndexS(:,2)], 'VariableNames',{'Y','Cond','Sub','Elec1','Elec2'});
                lmeTBL.Cond = nominal(lmeTBL.Cond);
                lmeTBL.Sub = nominal(lmeTBL.Sub);
                lmeTBL.Elec1 = nominal(lmeTBL.Elec1);
                lmeTBL.Elec2 = nominal(lmeTBL.Elec2);
                lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec1:Elec2)','fitmethod','reml','DummyVarCoding','effects');
                [~,~,lmeStats] = fixedEffects(lmeStruct);
                tMap(1,itime) = lmeStats.tStat(2);
                pMap(1,itime) = lmeStats.pValue(2);
                [randBeta,~,~] = randomEffects(lmeStruct);
                Z = designMatrix(lmeStruct,'random');
                Ycorr = lmeTBL.Y-Z*randBeta;
                obsVal1 = Ycorr(lmeTBL.Cond=='1');
                obsVal2 = Ycorr(lmeTBL.Cond=='2');
                y2plot(:,itime) = [mean(obsVal1);mean(obsVal2)];
                y2raw(:,itime) = [mean(lmeTBL.Y(lmeTBL.Cond=='1')),mean(lmeTBL.Y(lmeTBL.Cond=='2'))];
                se2plot(:,itime) = [std(obsVal1)./sqrt(numel(obsVal1));std(obsVal2)./sqrt(numel(obsVal2))];
                
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz' filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','y2raw','Para');
            
            
        end
        
        
        
        %% section4: calculate MVGC (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'mvgc')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(1,size(metricM,2));
            pMap = zeros(1,size(metricM,2));
            y2plot = zeros(2,size(metricM,2));
            se2plot = zeros(2,size(metricM,2));
            
            strlen = 0;
            % TFR point wise LME
            for itime = 1:size(metricM,2)
                
                
                s = ['Calculating freq point: ' num2str(itime) '/' num2str(size(metricM,2)) 'freq'];
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
                
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','Para');
            
            
        end
        
        
        %% section4-2: calculate Granger causality (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'Granger')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(1,size(allMetricM,2));
            pMap = zeros(1,size(allMetricM,2));
            y2plot = zeros(2,size(allMetricM,2));
            y2raw = zeros(2,size(allMetricM,2));
            se2plot = zeros(2,size(allMetricM,2));
            
            strlen = 0;
            % TFR point wise LME
            for ifreq = 1:size(allMetricM,2)
                
                
                s = ['Calculating freq point: ' num2str(ifreq) '/' num2str(size(allMetricM,2)) 'freqs'];
                strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
                strlen = strlentmp - strlen;
                
                frameDataM = double(squeeze(allMetricM(:,ifreq)));
                % skip nan point
                if ~any(frameDataM,'all')
                    continue
                end
                frameDataS = double(squeeze(allMetricS(:,ifreq)));
                
                
                %                 lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                %                     [elecIndexM;elecIndexS], 'VariableNames',{'Y','Cond','Sub','Elec'});
                %                 lmeTBL.Cond = nominal(lmeTBL.Cond);
                %                 lmeTBL.Sub = nominal(lmeTBL.Sub);
                %                 lmeTBL.Elec = nominal(lmeTBL.Elec);
                %                 lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
                
                lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                    [elecIndexM(:,1);elecIndexS(:,1)], [elecIndexM(:,2);elecIndexS(:,2)], 'VariableNames',{'Y','Cond','Sub','Elec1','Elec2'});
                lmeTBL.Cond = nominal(lmeTBL.Cond);
                lmeTBL.Sub = nominal(lmeTBL.Sub);
                lmeTBL.Elec1 = nominal(lmeTBL.Elec1);
                lmeTBL.Elec2 = nominal(lmeTBL.Elec2);
                lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec1:Elec2)','fitmethod','reml','DummyVarCoding','effects');
                
                [~,~,lmeStats] = fixedEffects(lmeStruct);
                tMap(1,ifreq) = lmeStats.tStat(2);
                pMap(1,ifreq) = lmeStats.pValue(2);
                [randBeta,~,~] = randomEffects(lmeStruct);
                Z = designMatrix(lmeStruct,'random');
                Ycorr = lmeTBL.Y-Z*randBeta;
                obsVal1 = Ycorr(lmeTBL.Cond=='1');
                obsVal2 = Ycorr(lmeTBL.Cond=='2');
                y2plot(:,ifreq) = [mean(obsVal1);mean(obsVal2)];
                y2raw(:,ifreq) = [mean(lmeTBL.Y(lmeTBL.Cond=='1')),mean(lmeTBL.Y(lmeTBL.Cond=='2'))];
                se2plot(:,ifreq) = [std(obsVal1)./sqrt(numel(obsVal1));std(obsVal2)./sqrt(numel(obsVal2))];
                
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','y2raw','Para');
            
            
            
            
        end
        
        
    end
end

varargout{1} = toc/3600;
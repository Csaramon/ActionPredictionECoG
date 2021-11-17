function varargout = ConnectivityCalculation(calculate)

tic
if nargin < 1
    calculate = 'COHtf'
end

% initialize base path and toolbox
if strcmpi(computer,'PCWIN64')
    addpath('C:\Users\qin2\Documents\MATLAB\toolbox')
    addpath('C:\Users\qin2\Documents\MATLAB\toolbox\fieldtrip-20210418')
    basePath = 'C:\Users\qin2\Documents\ActionPredictionECoG\';
elseif strcmpi(computer,'MACI64')
    addpath('~/Desktop/ActionPrediction')
    addpath('C:\Users\qin2\Documents\MATLAB\toolbox\fieldtrip-20210418')
    basePath = '~/Desktop/ActionPrediction/';
elseif strcmpi(computer,'GLNXA64')
    addpath('/data00/Chaoyi/ActionPredictionECoG/')
    addpath('/data00/Chaoyi/toolbox/fieldtrip-20210418/')
    basePath = '/data00/Chaoyi/ActionPredictionECoG/';
end

resultPath = [basePath 'Results' filesep];

ft_defaults

allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub



%% -------- ROI Partition --------

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
icontrol = [3];
allPair = nchoosek(seedIndex,2);
for iseed = seedIndex
    for isearch = searchIndex
        
        %                 skip redundant pairs
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
        
        
        % load altlas infomation for Control roi
        aparc_nii = load_nifti([basePath 'Atlas' filesep ROIAtlas{icontrol}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        if ndims(aparc_vol) == 3
            for iroi = ROIIndex{icontrol}
                temp_indx = find(aparc_vol==iroi);
                allroi_index = [allroi_index;temp_indx];
            end
        elseif ndims(aparc_vol) == 4
            %%%------ For probabilistic map ------%%%
            aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{icontrol}),4));
            allroi_index = find(aparc_vol>=10);  % minimum probability
        end
        
        % atlas parcellation coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        control_coordiantes = ras_coordiantes;
        
        
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
            subPath = [basePath filesep 'Data' filesep subname filesep];
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
            
            
            %% section2-3: calculate Coherence across time frequency (COHtf) %%
            %%%%%%%%%%%%%%%%%%%%%a%%%
            if strcmp(calculate,'COHtf')
                
                p=0.05; % threshold for IVC
                timeWin = 1; % unit in second
                timeStep = 0.1; % unit in second
                
                
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
                cfg.trials = find(trlData.trialinfo(:,1)==0);
                trlDataM = ft_selectdata(cfg,trlData);
                
                cfg = [];
                cfg.trials = find(trlData.trialinfo(:,1)==1);
                trlDataS = ft_selectdata(cfg,trlData);
                
                timeRange = [min(trlData.time{1}) max(trlData.time{1})];
                
                clear trlData
                
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
                    
                    % calculate fourier spectrum
                    cfg            = [];
                    cfg.output     = 'fourier';
                    cfg.method     = 'mtmfft';
                    cfg.foilim     = [2 120];
                    % cfg.foi          = logspace(log10(2),log10(128),32);
                    cfg.tapsmofrq  = 5;
                    cfg.keeptrials = 'yes';
                    freqM    = ft_freqanalysis(cfg, trlDataMtmp);

                    freqS    = ft_freqanalysis(cfg, trlDataStmp);
                    
                    % calculate COH in Matched Condition
                    cfg            = [];
                    cfg.method     = 'coh';
                    %                 cfg.complex = 'absimag';
                    cfg.channelcmb = {trlDataM.label(seedElec) trlDataM.label(searchElec)};
                    COHM             = ft_connectivityanalysis(cfg, freqM);
                    allcohM(:,:,in) = COHM.cohspctrm;
                    % calculate COH in Scambled Condition
                    cfg            = [];
                    cfg.method     = 'coh';
                    %                 cfg.complex = 'absimag';
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
            
            %% section3: calculate Phase Slope Index (PSI)%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'PSI')
                
                timeTFR = 0.002; % time step for time frequency results
                timeWin = 1; % unit in second
                timeStep = 0.05; % unit in second
                freqRange = [20 30];
                
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
                if ~any(seedElec) | ~any(searchElec) | isempty(setdiff(seedElec,searchElec))
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
                
                
                %%%%---------- hilbert ----------%%%%
                %                 cfg              = [];
                %                 cfg.output       = 'fourier';
                %                 cfg.method       = 'hilbert';
                %                 cfg.foi          = freqRange(1):5:freqRange(2);
                %                 cfg.width        =  5;
                %                 cfg.toi          = min(trlDataM.time{1}):timeTFR:max(trlDataM.time{1}); %'all';
                %                 cfg.keeptrials   = 'yes';
                %                 cfg.filttype = 'firws';
                %                 cfg.filtorder = nan;
                %                 cfg.filtdir = 'onepass-zerophase';
                %                 cfg.precision = 'single';
                
                %%%%---------- wavelet ----------%%%%
                cfg              = [];
                cfg.output       = 'fourier';
                cfg.method       = 'wavelet';
                cfg.foi          = freqRange(1):1:freqRange(2);
                cfg.width        =  5;
                cfg.toi          = 'all';
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
                
                % initialize MVGC toolbox
                try
                    run /data00/Chaoyi/toolbox/tools/mvgc_v1.0/startup.m
                catch
                    run C:\Users\qin2\Documents\MATLAB\toolbox\tools\mvgc_v1.0\startup.m
                end
                % time window use to calculate
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
                % choose control electrodes according to MNI coordinates
                tempdev = pdist2(elecposMNI,control_coordiantes);
                [it,~] = find(tempdev <=roiDist);
                controlElec = unique(it);
                
                % skip if no electrode pair
                if ~any(seedElec) | ~any(searchElec) | isempty(setdiff(seedElec,searchElec))
                    continue
                end
                
                cfg = [];
                cfg.demean = 'yes';
                cfg.detrend = 'yes';
                
                trlData = ft_preprocessing(cfg,trlData);
                % Downsample Data
                cfg = [];
                cfg.resamplefs = 250;
                cfg.detrend         = 'no';
                [trlData] = ft_resampledata(cfg, trlData);
                
                badChanInd = trlData.trial{1,1}(seedElec,1)==0;
                seedElec(badChanInd) = [];
                
                badChanInd = trlData.trial{1,1}(searchElec,1)==0;
                searchElec(badChanInd) = [];
                
                badChanInd = trlData.trial{1,1}(controlElec,1)==0;
                controlElec(badChanInd) = [];
                
                % remove superimposed electrodes
                duplicateInd = intersect(searchElec,seedElec);
                if ~isempty(duplicateInd)
                    if numel(searchElec) > numel(seedElec)
                        searchElec = setdiff(searchElec,duplicateInd);
                    else
                        seedElec = setdiff(seedElec,duplicateInd);
                    end
                end
                duplicateInd = intersect([searchElec;seedElec],controlElec);
                if ~isempty(duplicateInd)
                    controlElec = setdiff(controlElec,duplicateInd);
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
                clear trlData
                
                
                % calculate mvgc
                allGrangerM = [];
                allGrangerS = [];
                Npair = 1;
                allChanCmb = [];
                
                time2cal  = trlDataM.time{1}>=min(timeWin) & trlDataM.time{1}<=max(timeWin);
                
                controlElec = []; % to do pairwise
                if isempty(controlElec)
                    % calculate connectivity metric with two regions
                    for iseedElec = seedElec'
                        for isearchElec = searchElec'
                            allChanCmb = [allChanCmb;[iseedElec,isearchElec]];
                            % Matched Condition
                            X = [];
                            for itrl = 1:numel(trlDataM.trial)
                                X(:,:,itrl) = trlDataM.trial{itrl}([iseedElec,isearchElec],time2cal);
                            end
                            
                            mvgc_calculate % calculate granger causality
                            FM = f;
                            allGrangerM(Npair,:) = squeeze((FM(2,1,:)-FM(1,2,:))./ ...
                                (FM(1,2,:)+FM(2,1,:)))';
                            
                            % Scrambled Condition
                            X = [];
                            for itrl = 1:numel(trlDataS.trial)
                                X(:,:,itrl) = trlDataS.trial{itrl}([iseedElec,isearchElec],time2cal);
                            end
                            
                            mvgc_calculate % calculate granger causality
                            FS = f;
                            fres = size(f,3)-1;
                            freqPoints = sfreqs(fres,fs)';
                            allGrangerS(Npair,:) = squeeze((FS(2,1,:)-FS(1,2,:))./ ...
                                (FS(1,2,:)+FS(2,1,:)))';
                            
                            % count for pairs of eletrodes
                            Npair = Npair + 1;
                        end
                    end
                    
                else
                    % calculate connectivity metric with three regions
                    for iseedElec = seedElec'
                        for isearchElec = searchElec'
                            for icontrolElec = controlElec'
                                allChanCmb = [allChanCmb;[iseedElec,isearchElec]];
                                % Matched Condition
                                X = [];
                                for itrl = 1:numel(trlDataM.trial)
                                    X(:,:,itrl) = trlDataM.trial{itrl}([iseedElec,isearchElec,icontrolElec],time2cal);
                                end
                                
                                mvgc_calculate % calculate granger causality
                                FM = f;
                                allGrangerM(Npair,:) = squeeze((FM(2,1,:)-FM(1,2,:))./ ...
                                    (FM(1,2,:)+FM(2,1,:)))';
                                
                                % Scrambled Condition
                                X = [];
                                for itrl = 1:numel(trlDataS.trial)
                                    X(:,:,itrl) = trlDataS.trial{itrl}([iseedElec,isearchElec,icontrolElec],time2cal);
                                end
                                
                                mvgc_calculate % calculate granger causality
                                FS = f;
                                fres = size(f,3)-1;
                                freqPoints = sfreqs(fres,fs)';
                                allGrangerS(Npair,:) = squeeze((FS(2,1,:)-FS(1,2,:))./ ...
                                    (FS(1,2,:)+FS(2,1,:)))';
                                
                                % count for pairs of eletrodes
                                Npair = Npair + 1;
                            end
                        end
                    end
                    
                end
                
                % generate index for Subject Electrode and Trial
                
                %                 elecIndexM =  cat(1,elecIndexM,isub*1000+allChanCmb);
                elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                %                 elecIndexS = cat(1,elecIndexS,isub*1000+allChanCmb);
                elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,allGrangerM);
                allMetricS = cat(1,allMetricS,allGrangerS);
                
                
                Para.chanCMB{isub} = allChanCmb;
                Para.timeWin = timeWin;
                Para.freq = freqPoints;
                nelec = nelec+Npair-1;
                
                
            end
            
            
            %% section4-2: calculate Granger causality (non-parametric) %%
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'GrangerTF')
                
                p=0.05; % threshold for IVC
                timeWin = 1; % unit in second
                timeStep = 0.1; % unit in second
                
                
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
                cfg.trials = find(trlData.trialinfo(:,1)==0);
                trlDataM = ft_selectdata(cfg,trlData);
                
                cfg = [];
                cfg.trials = find(trlData.trialinfo(:,1)==1);
                trlDataS = ft_selectdata(cfg,trlData);
                
                timeRange = [min(trlData.time{1}) max(trlData.time{1})];
                
                clear trlData
                
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
                    cfg.foi     = 2:2:120; % need equidistant frequency bins for granger method
                    cfg.tapsmofrq  = 5;
                    cfg.keeptrials = 'yes';
                    cfg.pad='nextpow2';
                    ft_warning off
                    freqM    = ft_freqanalysis(cfg, trlDataMtmp);

                    freqS    = ft_freqanalysis(cfg, trlDataStmp);
                    
                    % calculate Granger Index in Matched Condition
                    grangercfg = [];
                    grangercfg.method  = 'granger';
                    grangercfg.granger.conditional = 'no';
                    grangercfg.granger.sfmethod = 'bivariate';
                    grangercfg.channelcmb = {trlDataM.label(seedElec) trlDataM.label(searchElec)};
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
                    allGrangerM(:,:,in) = (grangerM.grangerspctrm(fromInd,:)-grangerM.grangerspctrm(toInd,:))./ ...
                        (grangerM.grangerspctrm(fromInd,:)+grangerM.grangerspctrm(toInd,:));
                    
                    % calculate Granger Index in Scambled Condition
                    grangercfg = [];
                    grangercfg.method  = 'granger';
                    grangercfg.granger.conditional = 'no';
                    grangercfg.granger.sfmethod = 'bivariate';
                    grangercfg.channelcmb = {trlDataS.label(seedElec) trlDataS.label(searchElec)};
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
                    allGrangerS(:,:,in) = (grangerS.grangerspctrm(fromInd,:)-grangerS.grangerspctrm(toInd,:))./ ...
                        (grangerS.grangerspctrm(fromInd,:)+grangerS.grangerspctrm(toInd,:));
                    
                    
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
            
            
            %% section4-3: calculate Granger causality %%
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'Granger')
                
                p=0.05; % threshold for IVC
                timeWin = [-0.5 0]; % unit in second
                
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
                        
                        %                         % calculate multivariate autoregressive model
                        %                         cfg         = [];
                        %                         cfg.order   = 20;
                        %                         cfg.toolbox = 'bsmart';
                        %                         cfg.channel = [iseedElec,isearchElec];
                        %                         mdata       = ft_mvaranalysis(cfg, trlDataM);
                        %                         % calculate spectral transfer function
                        %                         cfg        = [];
                        %                         cfg.method = 'mvar';
                        %                         mfreq      = ft_freqanalysis(cfg, mdata);
                        %                         % calculate granger index
                        %                         cfg            = [];
                        %                         cfg.method     = 'granger';
                        %                         grangerM             = ft_connectivityanalysis(cfg, mfreq);
                        %                         if strcmp(trlData.label(iseedElec),grangerM.label(1))
                        %                             allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(1,2,:)-grangerM.grangerspctrm(2,1,:))./ ...
                        %                                 (grangerM.grangerspctrm(1,2,:)+grangerM.grangerspctrm(2,1,:)))';
                        %
                        %                         elseif strcmp(trlData.label(iseedElec),grangerM.label(2))
                        %                             allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(2,1,:)-grangerM.grangerspctrm(1,2,:))./ ...
                        %                                 (grangerM.grangerspctrm(1,2,:)+grangerM.grangerspctrm(2,1,:)))';
                        %                         end
                        
                        cfg = [];
                        cfg.method = 'mtmfft';
                        cfg.output = 'fourier';
                        cfg.channel = [iseedElec,isearchElec];
                        cfg.tapsmofrq = 5;
                        freqdataM = ft_freqanalysis(cfg, trlDataM);
                        
                        grangercfg = [];
                        grangercfg.method  = 'granger';
                        grangercfg.granger.conditional = 'no';
                        grangercfg.granger.sfmethod = 'bivariate';
                        
                        grangerM      = ft_connectivityanalysis(grangercfg, freqdataM);
                        
                        if strncmp(trlData.label{iseedElec},grangerM.labelcmb(1),length(trlData.label{iseedElec}))
                            allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(1,:)-grangerM.grangerspctrm(2,:))./ ...
                                (grangerM.grangerspctrm(1,:)+grangerM.grangerspctrm(2,:)))';
                        elseif strncmp(trlData.label{iseedElec},grangerM.labelcmb(2),length(trlData.label{iseedElec}))
                            allGrangerM(Npair,:) = squeeze((grangerM.grangerspctrm(2,:)-grangerM.grangerspctrm(1,:))./ ...
                                (grangerM.grangerspctrm(1,:)+grangerM.grangerspctrm(2,:)))';
                        end
                        
                        
                        % Scrambled Condition
                        
                        %                         % calculate multivariate autoregressive model
                        %                         cfg         = [];
                        %                         cfg.order   = 20;
                        %                         cfg.toolbox = 'bsmart';
                        %                         cfg.channel = [iseedElec,isearchElec];
                        %                         mdata       = ft_mvaranalysis(cfg, trlDataS);
                        %                         % calculate spectral transfer function
                        %                         cfg        = [];
                        %                         cfg.method = 'mvar';
                        %                         mfreq      = ft_freqanalysis(cfg, mdata);
                        %                         % calculate granger index
                        %                         cfg            = [];
                        %                         cfg.method     = 'granger';
                        %                         grangerS             = ft_connectivityanalysis(cfg, mfreq);
                        %                         if strcmp(trlData.label(iseedElec),grangerS.label(1))
                        %                             allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(1,2,:)-grangerS.grangerspctrm(2,1,:))./ ...
                        %                                 (grangerS.grangerspctrm(1,2,:)+grangerS.grangerspctrm(2,1,:)))';
                        %
                        %                         elseif strcmp(trlData.label(iseedElec),grangerS.label(2))
                        %                             allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(2,1,:)-grangerS.grangerspctrm(1,2,:))./ ...
                        %                                 (grangerS.grangerspctrm(1,2,:)+grangerS.grangerspctrm(2,1,:)))';
                        %                         end
                        cfg = [];
                        cfg.method = 'mtmfft';
                        cfg.output = 'fourier';
                        cfg.channel = [iseedElec,isearchElec];
                        cfg.tapsmofrq = 5;
                        freqdataS = ft_freqanalysis(cfg, trlDataS);
                        
                        grangercfg = [];
                        grangercfg.method  = 'granger';
                        grangercfg.granger.conditional = 'no';
                        grangercfg.granger.sfmethod = 'bivariate';
                        
                        grangerS      = ft_connectivityanalysis(grangercfg, freqdataS);
                        
                        if strncmp(trlData.label{iseedElec},grangerS.labelcmb(1),length(trlData.label{iseedElec}))
                            allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(1,:)-grangerS.grangerspctrm(2,:))./ ...
                                (grangerS.grangerspctrm(1,:)+grangerS.grangerspctrm(2,:)))';
                        elseif strncmp(trlData.label{iseedElec},grangerS.labelcmb(2),length(trlData.label{iseedElec}))
                            allGrangerS(Npair,:) = squeeze((grangerS.grangerspctrm(2,:)-grangerS.grangerspctrm(1,:))./ ...
                                (grangerS.grangerspctrm(1,:)+grangerS.grangerspctrm(2,:)))';
                        end
                        
                        
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
            
            
            %% section5: calculate cross regional PAC %%
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(calculate,'PAC')
                
                % calculation parameters
                pacMethod = 'mvl'; % coh,plv,mlv,mi,pac
                timeWin = [0 0.5];
                
                %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
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
                
                
                % additional preprocessing
                %             cfg = [];
                %             cfg.demean = 'yes';
                %             cfg.detrend = 'yes';
                %
                %             trlData = ft_preprocessing(cfg,trlData);
                
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
                cfg.trials = find(trlData.trialinfo(:,1)==0);
                trlDataM = ft_selectdata(cfg,trlData);
                
                cfg = [];
                cfg.trials = find(trlData.trialinfo(:,1)==1);
                trlDataS = ft_selectdata(cfg,trlData);
                
                
                % calculate single electrode PAC without shuffle
                
                Npair = 1;
                allChanCmb = [];
                
                foi          = logspace(log10(2),log10(128),32);
                width        =  logspace(log10(3),log10(10),32); % adjustive cycles
                
                %%%%---- calculate PAC in Intact condition ---- %%%%
                % time-frequency decomposition (wavelet)
                cfg              = [];
                cfg.channel = seedElec;
                cfg.output       = 'fourier';
                cfg.method       = 'wavelet';
                cfg.foi          = foi(foi>=1 & foi<=30);
                cfg.width        =  width(foi>=1 & foi<30);
                cfg.toi          = min(timeWin):0.002:max(timeWin);%'all';
                cfg.precision = 'single';
                
                ft_warning off
                freqlow = ft_freqanalysis(cfg,trlDataM);
                
                % time-frequency decomposition (wavelet)
                cfg              = [];
                cfg.channel = searchElec;
                cfg.output       = 'fourier';
                cfg.method       = 'wavelet';
                cfg.foi          = foi(foi>=30 );
                cfg.width        =  width(foi>=30);
                cfg.toi          = min(timeWin):0.002:max(timeWin);%'all';
                cfg.precision = 'single';
                
                ft_warning off
                freqhigh = ft_freqanalysis(cfg,trlDataM);
                
                cfg = [];
                cfg.method = pacMethod;
                cfg.chanlow = freqlow.label;
                cfg.chanhigh = freqhigh.label;
                cfg.keeptrials = 'no';
                crossfreq = ft_crossfrequencyanalysis(cfg, freqlow, freqhigh);
                allPACM = crossfreq.crsspctrm;
                
                %%%%---- calculate PAC in Scrambled condition ---- %%%%
                % time-frequency decomposition (wavelet)
                cfg              = [];
                cfg.channel = seedElec;
                cfg.output       = 'fourier';
                cfg.method       = 'wavelet';
                cfg.foi          = foi(foi>=1 & foi<=30);
                cfg.width        =  width(foi>=1 & foi<30);
                cfg.toi          = min(timeWin):0.002:max(timeWin);%'all';
                cfg.precision = 'single';
                
                ft_warning off
                freqlow = ft_freqanalysis(cfg,trlDataS);
                
                % time-frequency decomposition (wavelet)
                cfg              = [];
                cfg.channel = searchElec;
                cfg.output       = 'fourier';
                cfg.method       = 'wavelet';
                cfg.foi          = foi(foi>=30 );
                cfg.width        =  width(foi>=30);
                cfg.toi          = min(timeWin):0.002:max(timeWin);%'all';
                cfg.precision = 'single';
                
                ft_warning off
                freqhigh = ft_freqanalysis(cfg,trlDataS);
                
                cfg = [];
                cfg.method = pacMethod;
                cfg.chanlow = freqlow.label;
                cfg.chanhigh = freqhigh.label;
                cfg.keeptrials = 'no';
                crossfreq = ft_crossfrequencyanalysis(cfg, freqlow, freqhigh);
                allPACS = crossfreq.crsspctrm;
                
                allChanCmb = [allChanCmb;crossfreq.labelcmb];
                Npair = Npair + size(crossfreq.labelcmb,1);
                
                
                % generate index for Subject Electrode and Trial
                metricM = allPACM;
                metricS = allPACS;
                
                elecIndexM =  cat(1,elecIndexM,[nelec:nelec+Npair-2]');
                subIndexM = cat(1,subIndexM,repmat(isub,Npair-1,1));
                
                
                elecIndexS = cat(1,elecIndexS,[nelec:nelec+Npair-2]');
                subIndexS = cat(1,subIndexS,repmat(isub,Npair-1,1));
                
                % concontenate all trial responses in chosen ROI
                allMetricM = cat(1,allMetricM,metricM);
                allMetricS = cat(1,allMetricS,metricS);
                
                Para.chanCMB{isub} = allChanCmb;
                Para.pacMethod = pacMethod;
                Para.timeWin = timeWin;
                Para.freqlow = freqlow.freq;
                Para.freqhigh = freqhigh.freq;
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
            yraw = zeros(2,size(metricM,3));
            seraw = zeros(2,size(metricM,3));
            
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
                rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
                rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
                yraw(:,itime) = [mean(rawVal1);mean(rawVal2)];
                seraw(:,itime) = [std(rawVal1)./sqrt(numel(rawVal1));std(rawVal2)./sqrt(numel(rawVal2))];
                
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
        if strcmp(calculate,'COHtrl') | strcmp(calculate,'COHtf') | strcmp(calculate,'GrangerTF')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(size(metricM,2),size(metricM,3));
            pMap = zeros(size(metricM,2),size(metricM,3));
            y2plot = zeros(2,size(metricM,2),size(metricM,3));
            se2plot = zeros(2,size(metricM,2),size(metricM,3));
            yraw = zeros(2,size(metricM,2),size(metricM,3));
            seraw = zeros(2,size(metricM,2),size(metricM,3));
            
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
                    rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
                    rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
                    yraw(:,ifreq,itime) = [mean(rawVal1);mean(rawVal2)];
                    seraw(:,ifreq,itime) = [std(rawVal1)./sqrt(numel(rawVal1));std(rawVal2)./sqrt(numel(rawVal2))];
               
                end
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ROIText{iseed} '_'  ...
                ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','yraw','seraw','Para');
            
            
        end
        
        
        %% section3: calculate Phase Slope Index (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PSI')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(1,size(allMetricM,2));
            pMap = zeros(1,size(allMetricM,2));
            y2plot = zeros(2,size(allMetricM,2));
            se2plot = zeros(2,size(allMetricM,2));
            yraw = zeros(2,size(allMetricM,2));
            seraw = zeros(2,size(allMetricM,2));
            
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
                se2plot(:,itime) = [std(obsVal1)./sqrt(numel(obsVal1));std(obsVal2)./sqrt(numel(obsVal2))];
                rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
                rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
                yraw(:,itime) = [mean(rawVal1);mean(rawVal2)];
                seraw(:,itime) = [std(rawVal1)./sqrt(numel(rawVal1));std(rawVal2)./sqrt(numel(rawVal2))];
                
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz' filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','yraw','seraw','Para');
            
            
        end
        
        
        
        %% section4: calculate Granger causality (LMEM) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'Granger') | strcmp(calculate,'mvgc')
            
            if ~any(allMetricM) | ~any(allMetricS)
                continue
            end
            
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
                
                
                lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                    [elecIndexM;elecIndexS], 'VariableNames',{'Y','Cond','Sub','Elec'});
                lmeTBL.Cond = nominal(lmeTBL.Cond);
                lmeTBL.Sub = nominal(lmeTBL.Sub);
                lmeTBL.Elec = nominal(lmeTBL.Elec);
                lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
                
                %                                 lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                %                                     [elecIndexM(:,1);elecIndexS(:,1)], [elecIndexM(:,2);elecIndexS(:,2)], 'VariableNames',{'Y','Cond','Sub','Elec1','Elec2'});
                %                                 lmeTBL.Cond = nominal(lmeTBL.Cond);
                %                                 lmeTBL.Sub = nominal(lmeTBL.Sub);
                %                                 lmeTBL.Elec1 = nominal(lmeTBL.Elec1);
                %                                 lmeTBL.Elec2 = nominal(lmeTBL.Elec2);
                %                                 lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec1:Elec2)','fitmethod','reml','DummyVarCoding','effects');
                %
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
                rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
                rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
                yraw(:,ifreq) = [mean(rawVal1);mean(rawVal2)];
                seraw(:,ifreq) = [std(rawVal1)./sqrt(numel(rawVal1));std(rawVal2)./sqrt(numel(rawVal2))];
                
            end
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot', ...
                'yraw','seraw','Para');
            
            
            
            
        end
        
        %% section5: calculate cross regional PAC %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PAC')
            CondIndexM = ones(size(subIndexM));
            CondIndexS = 2*ones(size(subIndexS));
            
            tMap = zeros(size(metricM,2),size(metricM,3));
            pMap = zeros(size(metricM,2),size(metricM,3));
            y2plot = zeros(2,size(metricM,2),size(metricM,3));
            se2plot = zeros(2,size(metricM,2),size(metricM,3));
            yraw = zeros(2,size(metricM,2),size(metricM,3));
            seraw = zeros(2,size(metricM,2),size(metricM,3));
            
            strlen = 0;
            % point wise LME
            for ifreq = 1:size(metricM,2)
                for itime = 1:size(metricM,3)
                    
                    s = ['Calculating tf point: ' num2str(ifreq) '/' num2str(size(metricM,2)) 'AmpFreqs, ' ...
                        num2str(itime) '/' num2str(size(metricM,3)) 'PhFreqs'];
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
            
            if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4)])
            end
            save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) filesep ...
                ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot', ...
                'yraw','seraw','Para');
            
        end
        
        
        %             %% Save results
        %                         if ~exist([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'],'file')
        %                     mkdir([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'])
        %                 end
        %                 save([resultPath calculate filesep ROIAtlas{iseed}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz' ...
        %                     filesep ROIText{iseed} '_'  ROIText{isearch}],'lmeTBL','tMap','pMap','y2plot','se2plot','yraw','seraw','Para');
        
        
        
    end
end

varargout{1} = toc;

end

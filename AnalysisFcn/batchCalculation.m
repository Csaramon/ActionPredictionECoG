function varargout = batchCalculation(calculate)

tic;
if nargin < 1
    calculate = 'PACold'
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

% Anatomy20 ROIS
% ROIIndex = {[198],[209],[128],[181],[135]};
% ROIAtlas = {'Anatomy_v20.nii','Anatomy_v20.nii','Anatomy_v20.nii','Anatomy_v20.nii','Anatomy_v20.nii'};
% ROIText = {'V1','V2','V3','V4','V5'};


% ROIIndex = {[1008,2008],[1029,2029],[1031,2031],[1022,2022],[1024,2024], ...
%     [1028,2028],[1005,2005],[1011,2011]};
% ROIAtlas = {'aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii', ...
%     'aparc+aseg.nii','aparc+aseg.nii','aparc+aseg.nii'};
% ROIText = {'InferiorParietalLobe','SuperiorParietalLobe','SupraMarginal','Postcentral','Precentral', ...
%     'SuperiorFrontal','Cuneus','LateralOccipital'};
roiDist = 1; % maximum distance between electrodes and ROI voxels


for iatlas = [1,3,7]%[1,3,7,8]%1:numel(ROIIndex) %[1,3,7,8,9]
    
    % load altlas infomation
    aparc_nii = load_nifti([basePath 'Atlas' filesep ROIAtlas{iatlas}]);
    aparc_vol = round(aparc_nii.vol);
    allroi_index = [];
    %%%------ For maxprobabilistic map ------%%%
    if ndims(aparc_vol) == 3
        for iroi = ROIIndex{iatlas}
            temp_indx = find(aparc_vol==iroi);
            allroi_index = [allroi_index;temp_indx];
        end
    elseif ndims(aparc_vol) == 4
        %%%------ For probabilistic map ------%%%
        aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{iatlas}),4));
        allroi_index = find(aparc_vol>=10);  % minimum probability
    end
    
    % atlas parcellation coordinates
    [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
    ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
    ras_coordiantes = ras_coordiantes(1:3,:)';
    aparc_coordiantes = ras_coordiantes;
    
    
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
        
        %% section1: calculate ERP %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'ERP')
            
            % calculation parameters
            IVCselect = 0;pIVC=0.05; % Whether to use IVC selection; threshold for IVC
            
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LAR_trlData.mat'],'file')
                a = load([dataPath subname 'LAR_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = trlData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,jt] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            % implement electrode selection accroding to IVC
            if IVCselect
                if exist([dataPath subname 'IVC.mat'],'file')
                    load([dataPath subname 'IVC']);
                end
                % IVC selection criterion
                respElecInd = find(IVC.Intact.theta(:,2)<pIVC | IVC.Intact.alpha(:,2)<pIVC | IVC.Intact.beta(:,2)<pIVC | ...
                    IVC.Intact.lgamma(:,2)<pIVC | IVC.Intact.hgamma(:,2)<pIVC | IVC.Scamble.theta(:,2)<pIVC | ...
                    IVC.Scamble.alpha(:,2)<pIVC | IVC.Scamble.beta(:,2)<pIVC | IVC.Scamble.lgamma(:,2)<pIVC | IVC.Scamble.hgamma(:,2)<pIVC);
                
                ROIelec = intersect(ROIelec,respElecInd);% select electrodes indexes
            end
            
            
            cfg = [];
            cfg.demean = 'yes';
            trlData = ft_preprocessing(cfg,trlData);
            % seperate conditions
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==0);
            trlDataM = ft_selectdata(cfg,trlData);
            
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==1);
            trlDataS = ft_selectdata(cfg,trlData);
            
            % calculate timelock results (ERP)
            cfg = [];
            cfg.channel = 'all';
            cfg.keeptrials = 'yes';
            timelockM = ft_timelockanalysis(cfg,trlDataM);
            timelockS = ft_timelockanalysis(cfg,trlDataS);
            
            % generate index of Subject Electrode and Trial for LME
            metricM = timelockM.trial(:,ROIelec,:,:);
            metricS = timelockS.trial(:,ROIelec,:,:);
            
            elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricM,1))';
            elecIndexM = cat(1,elecIndexM,reshape(elecIndextmp,numel(elecIndextmp),1));
            subIndexM = cat(1,subIndexM,repmat(isub,numel(elecIndextmp),1));
            trlIndexM = cat(1,trlIndexM,repmat(ones(size(trlDataM.trialinfo(:,2))),size(metricM,2),1)); % uniform trials
            %             lmeInd = zeros(size(trlDataM.trialinfo(:,2)));
            %             for im = unique(trlDataM.trialinfo(:,2))'
            %                 a = find(trlDataM.trialinfo(:,2)==im);
            %                 rep = find(a(2:end)-a(1:end-1)>1);
            %                 if isempty(rep)
            %                     lmeInd(a) = a;
            %                 else
            %                     tmpInd = repmat(a(1:rep(1)),numel(rep)+1,1);
            %                     lmeInd(a) = tmpInd(1:numel(a));
            %                 end
            %             end
            %             trlIndexM = cat(1,trlIndexM,repmat(lmeInd,size(metricM,2),1));
            
            elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricS,1))';
            elecIndexS = cat(1,elecIndexS,reshape(elecIndextmp,numel(elecIndextmp),1));
            subIndexS = cat(1,subIndexS,repmat(isub,numel(elecIndextmp),1));
            trlIndexS = cat(1,trlIndexS,repmat(ones(size(trlDataS.trialinfo(:,2))),size(metricS,2),1)); % uniform trials
            %             lmeInd = zeros(size(trlDataS.trialinfo(:,2)));
            %             for im = unique(trlDataS.trialinfo(:,2))'
            %                 a = find(trlDataS.trialinfo(:,2)==im);
            %                 rep = find(a(2:end)-a(1:end-1)>1);
            %                 if isempty(rep)
            %                     lmeInd(a) = a;
            %                 else
            %                     tmpInd = repmat(a(1:rep(1)),numel(rep)+1,1);
            %                     lmeInd(a) = tmpInd(1:numel(a));
            %                 end
            %             end
            %
            %             trlIndexS = cat(1,trlIndexS,repmat(lmeInd+size(metricM,1),size(metricS,2),1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,reshape(metricM,prod([size(metricM,1),size(metricM,2)]),size(metricM,3),size(metricM,4)));
            allMetricS = cat(1,allMetricS,reshape(metricS,prod([size(metricS,1),size(metricS,2)]),size(metricS,3),size(metricS,4)));
            
            Para.elecposMNI = [Para.elecposMNI;timelockM.elec.elecposMNI(ROIelec,:)];
            Para.time = timelockM.time;
            Para.IVCselect = IVCselect;
            Para.pIVC = pIVC;
            nelec = nelec+numel(ROIelec);
            
            %%%%%%%%-------- Cluster based Permutation Test --------%%%%%%%%
            
            %             cfg = [];
            %             cfg.method            = 'montecarlo';           % use the Monte Carlo Method to calculate the significance probability
            %             cfg.statistic         = 'indepsamplesT';        % use the independent samples T-statistic as a measure to evaluate the effect at the sample level
            %             cfg.correctm          = 'cluster';
            %             cfg.clusteralpha      = clusterP;                   % alpha level of the sample-specific test statistic that will be used for thresholding
            %             cfg.clustertail       = 0;
            %             cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
            %             cfg.tail              = 0;                      % -1, 1 or 0 (default = 0); one-sided or two-sided test
            %             cfg.correcttail       = 'prob';                 % the two-sided test implies that we do non-parametric two tests
            %             cfg.alpha             = statP;                   % alpha level of the permutation test
            %             cfg.numrandomization  = 100;                   % number of draws from the permutation distribution
            %             cfg.design  = [1*ones(size(timelockM.trial,1),1);2*ones(size(timelockS.trial,1),1)];
            %             % cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
            %             cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
            %             cfg.channel           = 'all';
            %             cfg.neighbours        = [];
            
            %%%%%%%%-------- Cluster based Permutation Test --------%%%%%%%%
            
            
            %%%%%%%%-------- Independent T-test with Correction --------%%%%%%%%
            
            %         cfg = [];
            %         cfg.channel     = 'all'; %now all channels
            %         % cfg.latency     = [0.3 0.7];
            %         % cfg.avgovertime = 'yes';
            %         cfg.parameter   = 'trial';
            %         cfg.method      = 'analytic';
            %         cfg.statistic   = 'ft_statfun_indepsamplesT';
            %         cfg.alpha       = 0.05;
            %         cfg.correctm    = 'bonferroni';
            %
            %         cfg.design  = [1*ones(size(timelockM.trial,1),1);2*ones(size(timelockS.trial,1),1)];
            %         % cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
            %         cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
            %         % cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
            
            %%%%%%%%-------- Independent T-test with Correction --------%%%%%%%%
            
            
            %             timelockStat = ft_timelockstatistics(cfg, timelockM, timelockS);
            
        end
        
        
        %% section2: calculate ITC %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'ITC')
            
            statP = 0.05; % significance level for statistical test
            clusterP = 0.05; % significance level for cluster
            %             toi = -0.5:0.05:0.8;
            toi = 0:0.1:0.5;
            twin = 0.5;
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LAR_trlData.mat'],'file')
                a = load([dataPath subname 'LAR_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            % seperate conditions
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==0);
            trlDataM = ft_selectdata(cfg,trlData);
            
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==1);
            trlDataS = ft_selectdata(cfg,trlData);
            
            % datafolder setting
            datafolder = ['ITC' filesep subname filesep];
            if ~exist([resultPath datafolder],'file')
                mkdir(resultPath,datafolder)
            end
            
            tmpDatM = cell2mat(trlDataM.trial);
            trlDatM = reshape(tmpDatM,numel(trlDataM.label), numel(trlDataM.time{1}), []);
            
            tmpDatS = cell2mat(trlDataS.trial);
            trlDatS = reshape(tmpDatS,numel(trlDataS.label), numel(trlDataS.time{1}), []);
            
            CorrM = zeros(numel(trlDataM.label),numel(toi),nchoosek(numel(trlDataM.trial),2));
            CorrS = zeros(numel(trlDataS.label),numel(toi),nchoosek(numel(trlDataS.trial),2));
            CorrP = zeros(numel(trlDataS.label),numel(toi));
            maskP = zeros(numel(trlDataS.label),numel(toi));
            for ichan = 1:numel(trlData.label)
                
                chanDataM = squeeze(trlDatM(ichan,:,:));
                chanDataS = squeeze(trlDatS(ichan,:,:));
                
                for it = 1:numel(toi)
                    timeInd = trlData.time{1} >= toi(it) & trlData.time{1} <= (toi(it)+twin);
                    
                    chantimeDataM = chanDataM(timeInd,:);
                    if any(chantimeDataM(:))
                        chantimeCorrM = triu(corr(chantimeDataM),1);
                        CorrM(ichan,it,:) = chantimeCorrM(triu(ones(length(chantimeCorrM)),1)==1);
                        
                        chantimeDataS = chanDataS(timeInd,:);
                        chantimeCorrS= triu(corr(chantimeDataS),1);
                        CorrS(ichan,it,:) = chantimeCorrS(triu(ones(length(chantimeCorrS)),1)==1);
                        
                        [~,pp] = ttest2(CorrM(ichan,it,:) ,CorrS(ichan,it,:),'alpha',statP);
                        
                        CorrP(ichan,it) = pp;
                        maskP(ichan,it) = pp < (0.05/numel(toi));
                        
                    end
                end
                
                
            end
            
            cfg = [];
            cfg.channel = 'all';
            ITCM = ft_timelockanalysis(cfg,trlDataM);
            ITCM.time = toi+twin/2;
            ITCM.avg = squeeze(mean(CorrM,3));
            ITCM.var = nan(size(ITCM.avg));
            ITCM.dof = ones(size(ITCM.avg));
            ITCM.mask = maskP;
            
            ITCS = ft_timelockanalysis(cfg,trlDataS);
            ITCS.time = toi+twin/2;
            ITCS.avg = squeeze(mean(CorrS,3));
            ITCS.var = nan(size(ITCS.avg));
            ITCS.dof = ones(size(ITCS.avg));
            
            save([resultPath datafolder filesep 'ITCData'],'ITCM','ITCS','CorrP','twin');
            
            
        end
        
        
        
        %% section2-1: calculate ITC with permutation test%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'ITCperm')
            
            statP = 0.05; % significance level for statistical test
            shfTime = 1000;
            clusterP = 0.05; % significance level for cluster
            twin = [0 0.5];
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LAR_trlData.mat'],'file')
                a = load([dataPath subname 'LAR_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            
            % datafolder setting
            datafolder = ['ITCperm' filesep subname filesep];
            if ~exist([resultPath datafolder],'file')
                mkdir(resultPath,datafolder)
            end
            
            tmpDat = cell2mat(trlData.trial);
            trlDat = reshape(tmpDat,numel(trlData.label), numel(trlData.time{1}), []);
            CorrData = zeros(numel(trlData.label),nchoosek(numel(trlData.trial),2));
            CorrP = zeros(numel(trlData.label),1);
            CorrTstat = zeros(numel(trlData.label),1);
            shfStat = zeros(numel(trlData.label),shfTime);
            
            strlen = 0;
            
            for ichan = 1:numel(trlData.label)
                
                s = ['Calculating electrode: ' num2str(ichan) '/' num2str(numel(trlData.label))];
                strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
                strlen = strlentmp - strlen;
                
                chanData = squeeze(trlDat(ichan,:,:));
                timeInd = trlData.time{1} >= twin(1) & trlData.time{1} <= twin(2);
                
                chantimeData = chanData(timeInd,:);
                if any(chantimeData(:))
                    chantimeCorr = triu(corr(chantimeData),1);
                    tmpCorr = chantimeCorr(chantimeCorr~=0);
                    CorrData(ichan,:) = tmpCorr;
                    
                    [~,~,~,stat] = ttest(tmpCorr,0,'alpha',statP);
                    CorrTstat(ichan) = stat.tstat;
                    
                    
                    for ishf = 1:shfTime
                        data2shf = mat2cell(chantimeData',ones(size(chantimeData,2),1));
                        randindx =arrayfun(@(dummy) randi(size(chantimeData,1)), 1:size(chantimeData,2), 'UniformOutput', false)';
                        permData = cellfun(@(x,y) [x(y+1:end),x(1:y)], data2shf, randindx, 'UniformOutput',false);
                        permData = cell2mat(permData);
                        
                        shfCorr = triu(corr(permData'),1);
                        
                        [~,~,~,shfstat] = ttest(shfCorr(shfCorr~=0),0,'alpha',statP);
                        shfStat(ichan,ishf) = shfstat.tstat;
                        
                    end
                    
                    shfStatSort = sort(shfStat(ichan,:));
                    [~,minInd] = min(abs(stat.tstat-shfStatSort));
                    CorrP(ichan) = min([minInd/shfTime,1-minInd/shfTime]);
                end
                
                
            end
            
            
            
            save([resultPath datafolder filesep 'ITCpermData'],'CorrData','CorrP','CorrTstat','shfStat','twin');
            
            
        end
        
        
        
        %% section3: calculate power spectrum and peak frequency %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PeakFrequency')
            % set parameter
            freqRange = [1 128];
            timeWin = [0 1];
            timeStep = 0.01;
            
            % make data filefolder
            datafolder = ['PowerSpectrum' filesep '1-128Hz'];
            if ~exist([resultPath datafolder],'file')
                mkdir(resultPath,datafolder)
            end
            
            
            % load reref data
            if exist([dataPath subname 'LAR_trlData.mat'],'file')
                a = load([dataPath subname 'LAR_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            for ielec = 1:size(trlData.label,1)
                
                activecoordinates = trlData.elec.elecposMNI(ielec,:);
                
                cfg = [];
                cfg.channel           = trlData.label(ielec);
                [rerefSingle] = ft_selectdata(cfg, trlData);
                
                % time frequency analysis
                cfg              = [];
                cfg.output       = 'pow';
                cfg.method       = 'wavelet';
                cfg.foi          = logspace(log10(freqRange(1)),log10(freqRange(2)),128);
                cfg.width      = logspace(log10(3),log10(10),128); % adjustive cycles
                cfg.toi          = min(trlData.time{1}):timeStep:max(trlData.time{1}); %'all';
                cfg.pad = 'nextpow2';
                ft_warning off
                freq         = ft_freqanalysis(cfg,rerefSingle);
                
                % select data
                toi = cfg.toi>=timeWin(1) & cfg.toi<=timeWin(2);
                freq_pow = squeeze(freq.powspctrm(:,:,toi));
                
                if any(freq_pow(:))
                    spec_to_plot = nanmean(freq_pow,2)';
                    %             if strcmp(subname,'peijian') & ielec==79
                    %                 save([resultpath datafolder filesep subname '_elec' num2str(ielec)],'freq_b');
                    %             end
                    peak_frequency=getpeak(spec_to_plot,freq.freq,log10(freqRange(1)),log10(freqRange(2)),1);% find peak frequencies of each electrode
                    
                    
                    
                    ff = figure('Name',[subname '_powerspctrum_elec' num2str(ielec)]);
                    h = axes('Parent',ff);
                    hp = plot(h,log2(freq.freq),log2(spec_to_plot));
                    set(h,'xtick',[0:7],'xticklabel',[1,2,4,8,16,32,64,128]);
                    xlim([0 7])
                    
                    hold on;
                    for ipeak = 1:numel(peak_frequency)
                        temppow = log2(spec_to_plot(freq_b.freq==peak_frequency(ipeak)));
                        tempfreq = log2(peak_frequency(ipeak));
                        hl(ipeak) = plot(h,[tempfreq tempfreq],[0.95*temppow 1.05*temppow],'r');
                    end
                    
                    xlabel(h,'Frequency(log Hz)')
                    ylabel(h,'log Power(au.)')
                    title(h,'Power Spectrum')
                    
                    if any(peak_frequency)
                        setappdata(ff,'OSCstate',1);
                    else
                        setappdata(ff,'OSCstate',0);
                    end
                    setappdata(ff,'MNICoordinate',activecoordinates);
                    setappdata(ff,'power',spec_to_plot);
                    setappdata(ff,'frequency',freq.freq);
                    setappdata(ff,'peakFrequency',peak_frequency);
                    
                    
                    saveas(ff,[resultPath datafolder filesep subname '_elec' num2str(ielec)]);
                    close all
                end
                
            end
            
            
            
            
        end
        
        %% section3-2: calculate power spectrum difference %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PSD')
            
            % calculation parameters
            IVCselect = 0;pIVC=0.05; % Whether to use IVC selection; threshold for IVC
            
            timeWin = [0 1];
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LAR_trlData.mat'],'file')
                a = load([dataPath subname 'LAR_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = trlData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,jt] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            % implement electrode selection accroding to IVC
            if IVCselect
                if exist([dataPath subname 'IVC.mat'],'file')
                    load([dataPath subname 'IVC']);
                end
                % IVC selection criterion
                respElecInd = find(IVC.Intact.theta(:,2)<pIVC | IVC.Intact.alpha(:,2)<pIVC | IVC.Intact.beta(:,2)<pIVC | ...
                    IVC.Intact.lgamma(:,2)<pIVC | IVC.Intact.hgamma(:,2)<pIVC | IVC.Scamble.theta(:,2)<pIVC | ...
                    IVC.Scamble.alpha(:,2)<pIVC | IVC.Scamble.beta(:,2)<pIVC | IVC.Scamble.lgamma(:,2)<pIVC | IVC.Scamble.hgamma(:,2)<pIVC);
                
                ROIelec = intersect(ROIelec,respElecInd);% select electrodes indexes
            end
            
            
            cfg = [];
            cfg.demean = 'yes';
            cfg.detrend = 'yes';
            
            trlData = ft_preprocessing(cfg,trlData);
            
            % skip bad channels
            badChanInd = trlData.trial{1,1}(ROIelec,1)==0;
            ROIelec(badChanInd) = [];
            
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
            
            % calculate Power spectrum
            cfg            = [];
            cfg.output     = 'pow';
            cfg.method     = 'mtmfft';
            %             cfg.foilim     = [2 120];
            cfg.foi          = 2:1:120; %logspace(log10(2),log10(128),32);
            %                     cfg.taper      =  'hanning';
            cfg.tapsmofrq  = 5;
            cfg.keeptrials = 'no';
            cfg.pad='nextpow2';
            freqM    = ft_freqanalysis(cfg, trlDataM);
            
            freqS    = ft_freqanalysis(cfg, trlDataS);
            
            % generate index of Subject Electrode and Trial for LME
            metricM = freqM.powspctrm(ROIelec,:,:);
            metricS = freqS.powspctrm(ROIelec,:,:);
            
            %                 elecIndexM =  cat(1,elecIndexM,isub*1000+allChanCmb);
            elecIndexM =  cat(1,elecIndexM,[nelec:nelec+numel(ROIelec)-1]');
            subIndexM = cat(1,subIndexM,repmat(isub,numel(ROIelec),1));
            
            %                 elecIndexS = cat(1,elecIndexS,isub*1000+allChanCmb);
            elecIndexS = cat(1,elecIndexS,[nelec:nelec+numel(ROIelec)-1]');
            subIndexS = cat(1,subIndexS,repmat(isub,numel(ROIelec),1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,metricM);
            allMetricS = cat(1,allMetricS,metricS);
            
            Para.elecposMNI = [Para.elecposMNI;freqM.elec.elecposMNI(ROIelec,:)];
            Para.timeWin = timeWin;
            Para.freq = freqM.freq;
            
            nelec = nelec+numel(ROIelec);
        end
        
        
        %% section3-3: calculate power spectrum difference of all movie %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PSDmovie')
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'Bipolar_rerefData.mat'],'file')
                a = load([dataPath subname 'Bipolar_rerefData']);
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
            
            % downsampling data and eventdata
            cfg = [];
            cfg.resamplefs = 500;
            cfg.detrend         = 'yes';
            cfg.demean         = 'yes';
            rerefData = ft_resampledata(cfg, rerefData);
            fs = rerefData.fsample;
            camInfo(:,4) = num2cell(cell2mat(camInfo(:,4))./1000);
            
            % time-frequency decomposition (wavelet)
            %             cfg              = [];
            %             cfg.output       = 'pow';
            %             cfg.method       = 'wavelet';
            %             cfg.foi          = 2:1:120;
            %             cfg.width        =  7;
            %             cfg.toi          = min(rerefData.time{1}):0.1:max(rerefData.time{1});
            %             cfg.precision = 'single';
            %             cfg.pad='nextpow2';
            
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
            
            %                 elecIndexM =  cat(1,elecIndexM,isub*1000+allChanCmb);
            elecIndexM =  cat(1,elecIndexM,[nelec:nelec+numel(ROIelec)-1]');
            subIndexM = cat(1,subIndexM,repmat(isub,numel(ROIelec),1));
            
            %                 elecIndexS = cat(1,elecIndexS,isub*1000+allChanCmb);
            elecIndexS = cat(1,elecIndexS,[nelec:nelec+numel(ROIelec)-1]');
            subIndexS = cat(1,subIndexS,repmat(isub,numel(ROIelec),1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,metricM);
            allMetricS = cat(1,allMetricS,metricS);
            
            %---------- add trial level ----------%
            %             metricM = shiftdim(tsI,2);
            %             metricS = shiftdim(tsS,2);
            
            % generate index for Subject Electrode and Trial
            %             elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricM,1))';
            %             elecIndexM = cat(1,elecIndexM,reshape(elecIndextmp,numel(elecIndextmp),1));
            %             subIndexM = cat(1,subIndexM,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexM = cat(1,trlIndexM,repmat(ones(size(metricM,1),1),size(metricM,2),1)); % uniform trials
            %
            %             elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricS,1))';
            %             elecIndexS = cat(1,elecIndexS,reshape(elecIndextmp,numel(elecIndextmp),1));
            %             subIndexS = cat(1,subIndexS,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexS = cat(1,trlIndexS,repmat(ones(size(metricS,1),1),size(metricS,2),1)); % uniform trials
            %
            %             % concontenate all trial responses in chosen ROI
            %             allMetricM = cat(1,allMetricM,reshape(metricM,prod([size(metricM,1),size(metricM,2)]),size(metricM,3),size(metricM,4)));
            %             allMetricS = cat(1,allMetricS,reshape(metricS,prod([size(metricS,1),size(metricS,2)]),size(metricS,3),size(metricS,4)));
            %
            
            Para.elecposMNI = [Para.elecposMNI;freq.elec.elecposMNI(ROIelec,:)];
            Para.freq = cat(2,freq.freq,freq2.freq);
            
            nelec = nelec+numel(ROIelec);
        end
        
        %% section3-3: calculate power spectrum difference between views %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PSDviews')
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LAR_rerefData.mat'],'file')
                a = load([dataPath subname 'LAR_rerefData']);
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
            
            % downsampling data and eventdata
            cfg = [];
            cfg.resamplefs = 500;
            cfg.detrend         = 'yes';
            cfg.demean         = 'yes';
            rerefData = ft_resampledata(cfg, rerefData);
            fs = rerefData.fsample;
            camInfo(:,4) = num2cell(cell2mat(camInfo(:,4))./1000);
            
            % time-frequency decomposition (wavelet)
            %             cfg              = [];
            %             cfg.output       = 'pow';
            %             cfg.method       = 'wavelet';
            %             cfg.foi          = 2:1:120;
            %             cfg.width        =  7;
            %             cfg.toi          = min(rerefData.time{1}):0.1:max(rerefData.time{1});
            %             cfg.precision = 'single';
            %             cfg.pad='nextpow2';
            
            % time-frequency decomposition (multi-taper)
            cfg              = [];
            cfg.channel = ROIelec;
            cfg.output       = 'pow';
            cfg.method       = 'mtmconvol';
            cfg.foi          = 32:2:120;
            cfg.t_ftimwin  = 14./cfg.foi;
            cfg.tapsmofrq  = 0.3 *cfg.foi; % tapers=2*tw*fw-1
            cfg.toi          = min(rerefData.time{1}):0.1:max(rerefData.time{1});
            %                                 cfg.precision = 'single';
            cfg.pad='nextpow2';
            
            
            ft_warning off
            freq2 = ft_freqanalysis(cfg,rerefData);
            
            cfg.foi          = 2:1:30;
            cfg.t_ftimwin  = 4./cfg.foi;
            cfg.tapsmofrq  = 0.3 *cfg.foi; % tapers=2*tw*fw-1
            freq = ft_freqanalysis(cfg,rerefData);
            
            clear rerefData
            
            allTS = cat(2,freq.powspctrm,freq2.powspctrm);
            
            % seperate data view 1 and view 2
            if size(camInfo,1) >1590
                camInfo = camInfo(1:1590,:);
            end
            uniMovie = unique(cell2mat(camInfo(:,3)));
            view1Ind = [];
            view2Ind = [];
            for im = uniMovie'
                tempMovieInd = find(cell2mat(camInfo(:,3))==im);
                view1Ind = [view1Ind;tempMovieInd(1:numel(tempMovieInd)/2)];
                view2Ind = [view2Ind;tempMovieInd(numel(tempMovieInd)/2+1:end)];
            end
            
            if size(camInfo,1) <= 795
                camInfoView1 = camInfo;
                camInfoView2 = [];
            else
                camInfoView1 = camInfo(view1Ind,:);
                camInfoView2 = camInfo(view2Ind,:);
            end
            
            % predefine variable for the reconcantenate time series
            tsIview1 = [];
            tsSview1 = [];
            %%%%-------- reconcatenate data movie clip-wise and condition-wise --------%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % View 1
            for im = 1:size(camInfoView1,1)
                
                timeInd = freq.time >= camInfoView1{im,4} & freq.time <= (camInfoView1{im,4}+camInfoView1{im,6});
                % extract data of each movie
                tmpTS = allTS(:,:,timeInd);
                
                
                % concatenate data according to condition
                if camInfo{im,1}==0
                    tsIview1=cat(3,tsIview1,mean(tmpTS,3,'omitnan'));
                elseif camInfo{im,1}==1
                    tsSview1=cat(3,tsSview1,mean(tmpTS,3,'omitnan'));
                end
                
            end
            
            % View 2
            tsIview2 = [];
            tsSview2 = [];
            for im = 1:size(camInfoView2,1)
                
                timeInd = freq.time >= camInfoView2{im,4} & freq.time <= (camInfoView2{im,4}+camInfoView2{im,6});
                % extract data of each movie
                tmpTS = allTS(:,:,timeInd);
                
                
                % concatenate data according to condition
                if camInfo{im,1}==0
                    tsIview2=cat(3,tsIview2,mean(tmpTS,3,'omitnan'));
                elseif camInfo{im,1}==1
                    tsSview2=cat(3,tsSview2,mean(tmpTS,3,'omitnan'));
                end
                
            end
            
            % generate index of Subject Electrode and Trial for LME
            metricM = mean(tsSview1,3);
            metricS = mean(tsSview2,3);
            
            %                 elecIndexM =  cat(1,elecIndexM,isub*1000+allChanCmb);
            elecIndexM =  cat(1,elecIndexM,[nelec:nelec+size(metricM,1)-1]');
            subIndexM = cat(1,subIndexM,repmat(isub,size(metricM,1),1));
            
            %                 elecIndexS = cat(1,elecIndexS,isub*1000+allChanCmb);
            elecIndexS = cat(1,elecIndexS,[nelec:nelec+size(metricS,1)-1]');
            subIndexS = cat(1,subIndexS,repmat(isub,size(metricS,1),1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,metricM);
            allMetricS = cat(1,allMetricS,metricS);
            
            %---------- add trial level ----------%
            %             metricM = shiftdim(tsI,2);
            %             metricS = shiftdim(tsS,2);
            
            % generate index for Subject Electrode and Trial
            %             elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricM,1))';
            %             elecIndexM = cat(1,elecIndexM,reshape(elecIndextmp,numel(elecIndextmp),1));
            %             subIndexM = cat(1,subIndexM,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexM = cat(1,trlIndexM,repmat(ones(size(metricM,1),1),size(metricM,2),1)); % uniform trials
            %
            %             elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricS,1))';
            %             elecIndexS = cat(1,elecIndexS,reshape(elecIndextmp,numel(elecIndextmp),1));
            %             subIndexS = cat(1,subIndexS,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexS = cat(1,trlIndexS,repmat(ones(size(metricS,1),1),size(metricS,2),1)); % uniform trials
            %
            %             % concontenate all trial responses in chosen ROI
            %             allMetricM = cat(1,allMetricM,reshape(metricM,prod([size(metricM,1),size(metricM,2)]),size(metricM,3),size(metricM,4)));
            %             allMetricS = cat(1,allMetricS,reshape(metricS,prod([size(metricS,1),size(metricS,2)]),size(metricS,3),size(metricS,4)));
            %
            
            Para.elecposMNI = [Para.elecposMNI;freq.elec.elecposMNI(ROIelec,:)];
            Para.freq = cat(2,freq.freq,freq2.freq);
            
            nelec = nelec+numel(ROIelec);
        end
        
        
        %% section4: calculate time frequency representation %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'TFR')
            
            pIVC=0.05; % threshold for IVC
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            datafile = dir([resultPath calculate filesep subname filesep 'FreqData.mat']);
            load([datafile.folder filesep datafile.name]);
            
            if exist([dataPath subname 'IVC.mat'],'file')
                load([dataPath subname 'IVC']);
            end
            
            respElecInd = find(IVC.Intact.theta(:,2)<pIVC | IVC.Intact.alpha(:,2)<pIVC | IVC.Intact.beta(:,2)<pIVC | ...
                IVC.Intact.lgamma(:,2)<pIVC | IVC.Intact.hgamma(:,2)<pIVC | IVC.Scamble.theta(:,2)<pIVC | ...
                IVC.Scamble.alpha(:,2)<pIVC | IVC.Scamble.beta(:,2)<pIVC | IVC.Scamble.lgamma(:,2)<pIVC | IVC.Scamble.hgamma(:,2)<pIVC);
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = freqM.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,jt] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            ROIelec = intersect(ROIelec,respElecInd);
            
            
            % generate index for Subject Electrode and Trial
            %             metricM = freqM.powspctrm(:,ROIelec,:,:);
            %             metricS = freqS.powspctrm(:,ROIelec,:,:);
            metricM = abs(freqM.fourierspctrm(:,ROIelec,:,:)).^2;
            metricS = abs(freqS.fourierspctrm(:,ROIelec,:,:)).^2;
            
            elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricM,1))';
            elecIndexM = cat(1,elecIndexM,reshape(elecIndextmp,numel(elecIndextmp),1));
            subIndexM = cat(1,subIndexM,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexM =  cat(1,trlIndexM,repmat(freqM.trialinfo(:,2),size(metricM,2),1)); % matche movie
            %             trlIndexM = cat(1,trlIndexM,repmat([1:size(metricM,1)]',size(metricM,2),1)); % unique trials
            trlIndexM = cat(1,trlIndexM,repmat(ones(size(metricM,1),1),size(metricM,2),1)); % uniform trials
            
            %             lmeInd = zeros(size(freqM.trialinfo(:,2)));
            %             for im = unique(freqM.trialinfo(:,2))'
            %                 a = find(freqM.trialinfo(:,2)==im);
            %                 rep = find(a(2:end)-a(1:end-1)>1);
            %                 if isempty(rep)
            %                     lmeInd(a) = a;
            %                 else
            %                     tmpInd = repmat(a(1:rep(1)),numel(rep)+1,1);
            %                     lmeInd(a) = tmpInd(1:numel(a));
            %                 end
            %             end
            %             trlIndexM = cat(1,trlIndexM,repmat(lmeInd,size(metricM,2),1));
            
            elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricS,1))';
            elecIndexS = cat(1,elecIndexS,reshape(elecIndextmp,numel(elecIndextmp),1));
            subIndexS = cat(1,subIndexS,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexS = cat(1,trlIndexS,repmat(freqS.trialinfo(:,2),size(metricS,2),1)); % match movie
            %             trlIndexS = cat(1,trlIndexS,repmat([size(metricM,1)+1:size(metricM,1)+size(metricS,1)]',size(metricS,2),1)); % unique trials
            trlIndexS = cat(1,trlIndexS,repmat(ones(size(metricS,1),1),size(metricS,2),1)); % uniform trials
            
            %             lmeInd = zeros(size(freqS.trialinfo(:,2)));
            %             for im = unique(freqS.trialinfo(:,2))'
            %                 a = find(freqS.trialinfo(:,2)==im);
            %                 rep = find(a(2:end)-a(1:end-1)>1);
            %                 if isempty(rep)
            %                     lmeInd(a) = a;
            %                 else
            %                     tmpInd = repmat(a(1:rep(1)),numel(rep)+1,1);
            %                     lmeInd(a) = tmpInd(1:numel(a));
            %                 end
            %             end
            %
            %             trlIndexS = cat(1,trlIndexS,repmat(lmeInd+size(metricM,1),size(metricS,2),1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,reshape(metricM,prod([size(metricM,1),size(metricM,2)]),size(metricM,3),size(metricM,4)));
            allMetricS = cat(1,allMetricS,reshape(metricS,prod([size(metricS,1),size(metricS,2)]),size(metricS,3),size(metricS,4)));
            
            Para.elecposMNI = [Para.elecposMNI;freqM.elec.elecposMNI(ROIelec,:)];
            Para.time = freqM.time;
            Para.freq = freqM.freq;
            nelec = nelec+numel(ROIelec);
        end
        
        %% section4-2: calculate time frequency representation (elec level) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'TFRelec')
            
            pIVC=0.05; % threshold for IVC
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            datafile = dir([resultPath 'TFR' filesep subname filesep 'FreqData.mat']);
            load([datafile.folder filesep datafile.name]);
            
            if exist([dataPath subname 'IVC.mat'],'file')
                load([dataPath subname 'IVC']);
            end
            
            respElecInd = find(IVC.Intact.theta(:,2)<pIVC | IVC.Intact.alpha(:,2)<pIVC | IVC.Intact.beta(:,2)<pIVC | ...
                IVC.Intact.lgamma(:,2)<pIVC | IVC.Intact.hgamma(:,2)<pIVC | IVC.Scamble.theta(:,2)<pIVC | ...
                IVC.Scamble.alpha(:,2)<pIVC | IVC.Scamble.beta(:,2)<pIVC | IVC.Scamble.lgamma(:,2)<pIVC | IVC.Scamble.hgamma(:,2)<pIVC);
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = freqM.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,jt] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            %             ROIelec = intersect(ROIelec,respElecInd);
            
            
            % generate index for Subject Electrode and Trial
            %             metricM = freqM.powspctrm(:,ROIelec,:,:);
            %             metricS = freqS.powspctrm(:,ROIelec,:,:);
            metricM = squeeze(nanmean(abs(freqM.fourierspctrm(:,ROIelec,:,:)).^2,1));
            metricS = squeeze(nanmean(abs(freqS.fourierspctrm(:,ROIelec,:,:)).^2,1));
            
            elecIndexM =  cat(1,elecIndexM,[nelec:nelec+numel(ROIelec)-1]');
            subIndexM = cat(1,subIndexM,repmat(isub,numel(ROIelec),1));
            
            elecIndexS = cat(1,elecIndexS,[nelec:nelec+numel(ROIelec)-1]');
            subIndexS = cat(1,subIndexS,repmat(isub,numel(ROIelec),1));
            
            % concontenate all elec responses in chosen ROI
            allMetricM(nelec:nelec+numel(ROIelec)-1,:,:) = metricM;
            allMetricS(nelec:nelec+numel(ROIelec)-1,:,:) = metricS;
            
            
            Para.elecposMNI = [Para.elecposMNI;freqM.elec.elecposMNI(ROIelec,:)];
            Para.time = freqM.time;
            Para.freq = freqM.freq;
            nelec = nelec+numel(ROIelec);
        end
        
        
        %% section5: calculate band power %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'BP')
            
            timeWin = [-0.5 1];
            freqRange = [20 30] ; % [20 30] [60 90]
            
            %%%%%%%%%%%%%%% load trl data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            datafile = dir([dataPath subname 'LAR_trlData.mat']);
            load([datafile.folder filesep datafile.name]);
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = trlData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            cfg = [];
            cfg.demean = 'yes';
            cfg.detrend = 'yes';
            
            trlData = ft_preprocessing(cfg,trlData);
            
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
            
            % time-frequency decomposition (wavelet)
            %             cfg              = [];
            %             cfg.channel  = ROIelec;
            %             cfg.output       = 'pow';
            %             cfg.method       = 'wavelet';
            %             cfg.foi          = min(freqRange):max(freqRange);
            %             cfg.width        =  7;
            %             cfg.toi          = min(timeWin):0.01:max(timeWin);
            %             cfg.precision = 'single';
            %             cfg.keeptrials = 'yes';
            %             cfg.pad='nextpow2';
            
            % multi taper
            cfg              = [];
            cfg.channel = ROIelec;
            cfg.output       = 'pow';
            cfg.method       = 'mtmconvol';
            cfg.foi          = min(freqRange):max(freqRange);
            cfg.t_ftimwin  = 4./cfg.foi;
            cfg.tapsmofrq  = 0.3 *cfg.foi; % tapers=2*tw*fw-1
            cfg.toi          = 'all';%min(timeWin):0.01:max(timeWin);
            cfg.keeptrials = 'yes';
            %                                 cfg.precision = 'single';
            cfg.pad='nextpow2';
            
            ft_warning off
            freqM = ft_freqanalysis(cfg,trlDataM);
            
            freqS = ft_freqanalysis(cfg,trlDataS);
            
            % raw power
            
            % normalise power
            tmpPowerM = freqM.powspctrm;
            tmpPowerS =freqS.powspctrm;
            %             minPower = min(min(cat(1,tmpPowerM,tmpPowerS),[],1),[],4);
            %             maxPower = max(max(cat(1,tmpPowerM,tmpPowerS),[],1),[],4);
            %             minPowerM = repmat(minPower,size(tmpPowerM,1),1,1,size(tmpPowerM,4));
            %             maxPowerM = repmat(maxPower,size(tmpPowerM,1),1,1,size(tmpPowerM,4));
            %             minPowerS = repmat(minPower,size(tmpPowerS,1),1,1,size(tmpPowerS,4));
            %             maxPowerS = repmat(maxPower,size(tmpPowerS,1),1,1,size(tmpPowerS,4));
            
            %             metricM = (tmpPowerM-minPowerM)./(maxPowerM-minPowerM);
            metricM = tmpPowerM;
            metricM = nanmean(metricM,3);
            
            
            %             metricS = (tmpPowerS-minPowerS)./(maxPowerS-minPowerS);
            metricS = tmpPowerS;
            metricS = nanmean(metricS,3);
            
            
            %---------- elec level ----------%
            metricM = squeeze(mean(metricM,1));
            metricS = squeeze(mean(metricS,1));
            if size(metricM,1) > size(metricM,2)
                metricM = metricM';
                metricS =metricS';
            end
            %                 elecIndexM =  cat(1,elecIndexM,isub*1000+allChanCmb);
            elecIndexM =  cat(1,elecIndexM,[nelec:nelec+numel(ROIelec)-1]');
            subIndexM = cat(1,subIndexM,repmat(isub,numel(ROIelec),1));
            
            %                 elecIndexS = cat(1,elecIndexS,isub*1000+allChanCmb);
            elecIndexS = cat(1,elecIndexS,[nelec:nelec+numel(ROIelec)-1]');
            subIndexS = cat(1,subIndexS,repmat(isub,numel(ROIelec),1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,metricM);
            allMetricS = cat(1,allMetricS,metricS);
            
            %---------- add trial level ----------%
            % generate index for Subject Electrode and Trial
            %             elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricM,1))';
            %             elecIndexM = cat(1,elecIndexM,reshape(elecIndextmp,numel(elecIndextmp),1));
            %             subIndexM = cat(1,subIndexM,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexM = cat(1,trlIndexM,repmat(ones(size(metricM,1),1),size(metricM,2),1)); % uniform trials
            %
            %
            %             elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricS,1))';
            %             elecIndexS = cat(1,elecIndexS,reshape(elecIndextmp,numel(elecIndextmp),1));
            %             subIndexS = cat(1,subIndexS,repmat(isub,numel(elecIndextmp),1));
            %             trlIndexS = cat(1,trlIndexS,repmat(ones(size(metricS,1),1),size(metricS,2),1)); % uniform trials
            %
            %
            %             % concontenate all trial responses in chosen ROI
            %             allMetricM = cat(1,allMetricM,reshape(metricM,prod([size(metricM,1),size(metricM,2)]),size(metricM,3),size(metricM,4)));
            %             allMetricS = cat(1,allMetricS,reshape(metricS,prod([size(metricS,1),size(metricS,2)]),size(metricS,3),size(metricS,4)));
            
            
            Para.elecposMNI = [Para.elecposMNI;freqM.elec.elecposMNI(ROIelec,:)];
            Para.time = freqM.time;
            Para.freq = freqM.freq;
            nelec = nelec+numel(ROIelec);
            
            
        end
        
        %% section5-2: calculate band power of all movies%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'BPall')
            
            freqRange = [13 30];
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LAR_rerefData.mat'],'file')
                a = load([dataPath subname 'LAR_rerefData']);
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
            
            % downsampling data and eventdata
            cfg = [];
            cfg.resamplefs = 500;
            cfg.detrend         = 'yes';
            cfg.demean         = 'yes';
            rerefData = ft_resampledata(cfg, rerefData);
            fs = rerefData.fsample;
            camInfo(:,4) = num2cell(round(cell2mat(camInfo(:,4))./1000.*fs));
            
            foi          = logspace(log10(2),log10(128),32);
            width        =  logspace(log10(3),log10(10),32); % adjustive cycles
            
            % time-frequency decomposition (wavelet)
            cfg              = [];
            cfg.output       = 'pow';
            cfg.method       = 'wavelet';
            cfg.foi          = foi(foi>=freqRange(1) & foi<=freqRange(2));
            cfg.width        =  width(foi>=freqRange(1) & foi<=freqRange(2));
            cfg.toi          = 'all';
            cfg.precision = 'single';
            
            ft_warning off
            freq = ft_freqanalysis(cfg,rerefData);
            
            allTS = freq.powspctrm;
            % low pass filtering
            for ifreq = 1:size(allTS,2)
                [allTS(:,ifreq,:), ~, ~] = ft_preproc_lowpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
            end
            
            % calculate MPD
            MeanPowDiff
            
            metricM = tsI(ROIelec,:)';
            metricS = tsS(ROIelec,:)';
            
            % generate index for Subject Electrode and Trial
            elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricM,1))';
            elecIndexM = cat(1,elecIndexM,reshape(elecIndextmp,numel(elecIndextmp),1));
            subIndexM = cat(1,subIndexM,repmat(isub,numel(elecIndextmp),1));
            trlIndexM = cat(1,trlIndexM,repmat(ones(size(metricM,1),1),size(metricM,2),1)); % uniform trials
            
            
            elecIndextmp = repmat([nelec:nelec+numel(ROIelec)-1]',1,size(metricS,1))';
            elecIndexS = cat(1,elecIndexS,reshape(elecIndextmp,numel(elecIndextmp),1));
            subIndexS = cat(1,subIndexS,repmat(isub,numel(elecIndextmp),1));
            trlIndexS = cat(1,trlIndexS,repmat(ones(size(metricS,1),1),size(metricS,2),1)); % uniform trials
            
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,reshape(metricM,prod([size(metricM,1),size(metricM,2)]),size(metricM,3),size(metricM,4)));
            allMetricS = cat(1,allMetricS,reshape(metricS,prod([size(metricS,1),size(metricS,2)]),size(metricS,3),size(metricS,4)));
            
            Para.elecposMNI = [Para.elecposMNI;rerefData.elec.elecposMNI(ROIelec,:)];
            Para.freq = freqRange;
            nelec = nelec+numel(ROIelec);
            
            
        end
        
        
        %% section5-3: calculate electrode averaged power %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'avgpow')
            
            freqRange = [3 8];
            timeWin = [0 1];
            %%%%%%%%%%%%%%% load freq data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            datafile = dir([resultPath 'TFR' filesep subname filesep 'FreqData.mat']);
            load([datafile.folder filesep datafile.name]);
            
            
            % extract frequency points of intrest
            freqPoint = find(freqM.freq>=freqRange(1) & freqM.freq<=freqRange(2));
            timePoint = find(freqM.time>=timeWin(1) & freqM.time<=timeWin(2));
            % normalise power
            tmpPowerM = abs(freqM.fourierspctrm(:,:,freqPoint,:)).^2;
            tmpPowerS = abs(freqS.fourierspctrm(:,:,freqPoint,:)).^2;
            meanPower = (nanmean(nanmean(tmpPowerM,1),4)+nanmean(nanmean(tmpPowerS,1),4))/2;
            
            metricM = tmpPowerM./repmat(meanPower,size(tmpPowerM,1),1,1,size(tmpPowerM,4));
            metricM = nanmean(nanmean(metricM(:,:,:,timePoint),3),4);
            
            metricS = tmpPowerS./repmat(meanPower,size(tmpPowerS,1),1,1,size(tmpPowerS,4));
            metricS = nanmean(nanmean(metricS(:,:,:,timePoint),3),4);
            
            [~,pIvS] = ttest2(metricM,metricS);
            
            Para.time = timeWin;
            Para.freq = freqRange;
            
            if ~exist([resultPath calculate filesep ROIAtlas{1}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'],'file')
                mkdir([resultPath calculate filesep ROIAtlas{1}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz'])
            end
            save([resultPath calculate filesep ROIAtlas{1}(1:end-4) num2str(freqRange(1)) '_' num2str(freqRange(2)) 'Hz' filesep subname], ...
                'metricM','metricS','pIvS','Para');
            
            
        end
        
        
        %% section6: calculate PAC (bestPAC) %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PACold--')
            
            % calculation parameters
            timeWin = [0 1];
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LARER_trlData.mat'],'file')
                a = load([dataPath subname 'LARER_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = trlData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            % skip bad channels
            badChanInd = trlData.trial{1,1}(ROIelec,1)==0;
            ROIelec(badChanInd) = [];
            
            % skip if no electrode pair
            if ~any(ROIelec)
                continue
            end
            
            cfg = [];
            cfg.demean = 'yes';
            cfg.detrend = 'yes';
            
            trlData = ft_preprocessing(cfg,trlData);
            
            % seperate conditions
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==0);
            trlDataM = ft_selectdata(cfg,trlData);
            
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==1);
            trlDataS = ft_selectdata(cfg,trlData);
            
            
            cfg = [];
            cfg.channel = 'all';
            cfg.keeptrials = 'yes';
            timelockM = ft_timelockanalysis(cfg,trlDataM);
            timelockS = ft_timelockanalysis(cfg,trlDataS);
            
            if ~exist([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}])
            end
            
            
            % calculate single electrode PAC with shuffle
            timeInd = timelockM.time>=timeWin(1) & timelockM.time<=timeWin(2);
            allPACM = [];
            allPACS = [];
            Npair = 1;
            for ie = ROIelec'
                metricM = squeeze(timelockM.trial(:,ie,timeInd));
                metricS = squeeze(timelockS.trial(:,ie,timeInd));
                
                [pacmat, freqvec_ph, freqvec_amp] = find_pac_shf(metricM', trlDataM.fsample, 'mi');
                setappdata(gcf,'pacmat',pacmat);
                setappdata(gcf,'freqvec_ph',freqvec_ph);
                setappdata(gcf,'freqvec_amp',freqvec_amp);
                saveas(gcf,[resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas} filesep ...
                    subname num2str(ie) 'M.fig']);
                close(gcf)
                allPACM(Npair,:,:) = pacmat;
                
                [pacmat, freqvec_ph, freqvec_amp] = find_pac_shf(metricS', trlDataS.fsample, 'mi');
                setappdata(gcf,'pacmat',pacmat);
                setappdata(gcf,'freqvec_ph',freqvec_ph);
                setappdata(gcf,'freqvec_amp',freqvec_amp);
                saveas(gcf,[resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas} filesep ...
                    subname num2str(ie) 'S.fig']);
                close(gcf)
                allPACS(Npair,:,:) = pacmat;
                
                % count for pairs of eletrodes
                Npair = Npair + 1;
                
            end
            
            
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
            
            Para.chanCMB{isub} = ROIelec;
            Para.freqvec_ph = freqvec_ph;
            Para.freqvec_amp = freqvec_amp;
            nelec = nelec+Npair-1;
            
        end
        
        %% section6: calculate PAC (circular) %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PACold')
            
            % calculation parameters
            timeWin = [0 1];
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LARER_trlData.mat'],'file')
                a = load([dataPath subname 'LARER_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = trlData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            % skip bad channels
            badChanInd = trlData.trial{1,1}(ROIelec,1)==0;
            ROIelec(badChanInd) = [];
            
            % skip if no electrode pair
            if ~any(ROIelec)
                continue
            end
            
            cfg = [];
            cfg.demean = 'yes';
            cfg.detrend = 'yes';
            
            trlData = ft_preprocessing(cfg,trlData);
            
            % seperate conditions
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==0);
            trlDataM = ft_selectdata(cfg,trlData);
            
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==1);
            trlDataS = ft_selectdata(cfg,trlData);
            
            if ~exist([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}],'file')
                mkdir([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}])
            end
            
            % calculate Power spectrum
            cfg            = [];
            cfg.output     = 'fourier';
            cfg.method     = 'wavelet';
            cfg.channel = ROIelec;
            cfg.foi          = [2:2:30,35:5:120]; %logspace(log10(2),log10(128),32);
            cfg.toi = min(timeWin):0.002:max(timeWin);
            cfg.width = 5;
            cfg.keeptrials = 'yes';
            
            ft_warning off
            freqM    = ft_freqanalysis(cfg, trlDataM);
            freqS    = ft_freqanalysis(cfg, trlDataS);
            
            % calculate single electrode PAC with shuffle
            allPACM = [];
            allPACS = [];
            Npair = 1;
            for ie = 1:numel(ROIelec)
                
                [pacmat, pacmatSig] = find_pac(freqM, freqM, [ie,ie]);
                setappdata(gcf,'pacmat',pacmat);
                setappdata(gcf,'pacmatSig',pacmatSig);
                setappdata(gcf,'freqlow',freqM.freq(freqM.freq<=30));
                setappdata(gcf,'freqhigh',freqM.freq(freqM.freq>30));
                saveas(gcf,[resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas} filesep ...
                    subname num2str(ROIelec(ie)) 'M.fig']);
                close(gcf)
                allPACM(Npair,:,:) = pacmat;
                
                [pacmat, pacmatSig] = find_pac(freqS, freqS, [ie,ie]);
                setappdata(gcf,'pacmat',pacmat);
                setappdata(gcf,'pacmatSig',pacmatSig);
                setappdata(gcf,'freqlow',freqS.freq(freqS.freq<=30));
                setappdata(gcf,'freqhigh',freqS.freq(freqS.freq>30));
                saveas(gcf,[resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas} filesep ...
                    subname num2str(ROIelec(ie)) 'S.fig']);
                close(gcf)
                allPACS(Npair,:,:) = pacmat;
                
                % count for pairs of eletrodes
                Npair = Npair + 1;
                
            end
            
            
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
            
            Para.chanCMB{isub} = ROIelec;
            Para.freqlow = freqM.freq(freqM.freq<=30);
            Para.freqhigh = freqM.freq(freqM.freq>30);
            nelec = nelec+Npair-1;
            
        end
        
        
        %% section6: calculate PAC %%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(calculate,'PAC')
            
            % calculation parameters
            pacMethod = 'mvl'; % coh,plv,mlv,mi,pac
            timeWin = [0 1];
            
            %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist([dataPath subname 'LARER_trlData.mat'],'file')
                a = load([dataPath subname 'LARER_trlData']);
            end
            c = fieldnames(a);
            trlData = a.(c{1});
            
            % choose ROI electrodes according to MNI coordinates
            elecposMNI = trlData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,aparc_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            ROIelec = unique(it);
            
            % skip bad channels
            badChanInd = trlData.trial{1,1}(ROIelec,1)==0;
            ROIelec(badChanInd) = [];
            
            % skip if no electrode pair
            if ~any(ROIelec)
                continue
            end
            
            cfg = [];
            cfg.demean = 'yes';
            cfg.detrend = 'yes';
            
            trlData = ft_preprocessing(cfg,trlData);
            
            % seperate conditions
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==0);
            trlDataM = ft_selectdata(cfg,trlData);
            
            cfg = [];
            cfg.trials = find(trlData.trialinfo(:,1)==1);
            trlDataS = ft_selectdata(cfg,trlData);
            
            
            %%%%---- calculate PAC in Intact condition ---- %%%%
            % calculate Power spectrum
            cfg            = [];
            cfg.output     = 'fourier';
            cfg.method     = 'wavelet';
            cfg.channel = ROIelec;
            cfg.foi          = [2:2:30,35:5:120]; %logspace(log10(2),log10(128),32);
            cfg.toi = min(timeWin):0.002:max(timeWin);
            cfg.width = 5;
            cfg.keeptrials = 'yes';
            
            ft_warning off
            freqM    = ft_freqanalysis(cfg, trlDataM);
            
            freqS    = ft_freqanalysis(cfg, trlDataS);
            
            % calculate electrode PAC
            
            cfg = [];
            cfg.method = pacMethod;
            cfg.channel = 'all';
            cfg.freqlow = [2 30];
            cfg.freqhigh = [31 120];
            cfg.keeptrials = 'no';
            crossfreqM = ft_crossfrequencyanalysis(cfg, freqM);
            
            crossfreqS = ft_crossfrequencyanalysis(cfg, freqS);
            
            % generate index for Subject Electrode and Trial
            metricM = crossfreqM.crsspctrm;
            metricS = crossfreqS.crsspctrm;
            
            elecIndexM =  cat(1,elecIndexM,[nelec:nelec+numel(ROIelec)-1]');
            subIndexM = cat(1,subIndexM,repmat(isub,numel(ROIelec),1));
            
            elecIndexS = cat(1,elecIndexS,[nelec:nelec+numel(ROIelec)-1]');
            subIndexS = cat(1,subIndexS,repmat(isub,numel(ROIelec),1));
            
            % concontenate all trial responses in chosen ROI
            allMetricM = cat(1,allMetricM,metricM);
            allMetricS = cat(1,allMetricS,metricS);
            
            Para.chanCMB{isub} = ROIelec;
            Para.freqhigh = crossfreqM.freqhigh;
            Para.freqlow = crossfreqM.freqlow;
            nelec = nelec+numel(ROIelec);
            
        end
        
        
        
    end
    
    
    
    
    %% -------- Group Level --------
    
    %% section1: calculate event related potiential %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'ERP')
        CondIndexM = ones(size(subIndexM));
        CondIndexS = 2*ones(size(subIndexS));
        allMetricM = squeeze(allMetricM);
        allMetricS = squeeze(allMetricS);
        
        tMap = zeros(1,size(metricM,3));
        pMap = zeros(1,size(metricM,3));
        y2plot = zeros(2,size(metricM,3));
        se2plot = zeros(2,size(metricM,3));
        yraw = zeros(2,size(allMetricM,3));
        seraw = zeros(2,size(allMetricM,3));
        
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
        
        if ~exist([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)],'file')
            mkdir([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)])
        end
        save([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}], ...
            'lmeTBL','tMap','pMap','y2plot','se2plot','yraw','seraw','Para');
        
        
    end
    
    
    %% section3-2: calculate power spectrum difference %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'PSD') | strcmp(calculate,'PSDmovie') | strcmp(calculate,'PSDviews')
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
            
            
            %             lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
            %                 [elecIndexM;elecIndexS],[trlIndexM;trlIndexS], 'VariableNames',{'Y','Cond','Sub','Elec','Trl'});
            %             lmeTBL.Cond = nominal(lmeTBL.Cond);
            %             lmeTBL.Sub = nominal(lmeTBL.Sub);
            %             lmeTBL.Elec = nominal(lmeTBL.Elec);
            %             lmeTBL.Trl = nominal(lmeTBL.Trl);
            %             lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)+(1|Trl)','fitmethod','reml','DummyVarCoding','effects');
            
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
        
        if ~exist([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)],'file')
            mkdir([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)])
        end
        save([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}], ...
            'allMetricM','allMetricS','lmeTBL','tMap','pMap','y2plot','se2plot','yraw','seraw','Para');
        
        
    end
    
    %% section4: calculate time frequency representation %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'TFR')
        CondIndexM = ones(size(subIndexM));
        CondIndexS = 2*ones(size(subIndexS));
        
        tMap = zeros(size(metricM,3),size(metricM,4));
        pMap = zeros(size(metricM,3),size(metricM,4));
        y2plot = zeros(2,size(metricM,3),size(metricM,4));
        se2plot = zeros(2,size(metricM,3),size(metricM,4));
        
        strlen = 0;
        % TFR point wise LME
        for ifreq = 1:size(metricM,3)
            for itime = 1:size(metricM,4)
                
                s = ['Calculating tf point: ' num2str(ifreq) '/' num2str(size(metricM,3)) 'freqs, ' ...
                    num2str(itime) '/' num2str(size(metricM,4)) 'times'];
                strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
                strlen = strlentmp - strlen;
                
                frameDataM = double(squeeze(allMetricM(:,ifreq,itime)));
                % skip nan point
                if ~any(frameDataM,'all')
                    continue
                end
                frameDataS = double(squeeze(allMetricS(:,ifreq,itime)));
                
                
                lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
                    [elecIndexM;elecIndexS],[trlIndexM;trlIndexS], 'VariableNames',{'Y','Cond','Sub','Elec','Trl'});
                lmeTBL.Cond = nominal(lmeTBL.Cond);
                lmeTBL.Sub = nominal(lmeTBL.Sub);
                lmeTBL.Elec = nominal(lmeTBL.Elec);
                lmeTBL.Trl = nominal(lmeTBL.Trl);
                lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)+(1|Trl)','fitmethod','reml','DummyVarCoding','effects');
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
        
        if ~exist([resultPath calculate filesep ROIAtlas{1}(1:end-4)],'file')
            mkdir([resultPath calculate filesep ROIAtlas{1}(1:end-4)])
        end
        save([resultPath calculate filesep ROIAtlas{1}(1:end-4) filesep ROIText{iatlas}],'lmeTBL','tMap','pMap','y2plot','se2plot','Para');
        
        
    end
    
    %% section4-2: calculate time frequency representation (elec level) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'TFRelec')
        CondIndexM = ones(size(subIndexM));
        CondIndexS = 2*ones(size(subIndexS));
        
        tMap = zeros(size(allMetricM,2),size(allMetricM,3));
        pMap = zeros(size(allMetricM,2),size(allMetricM,3));
        y2plot = zeros(2,size(allMetricM,2),size(allMetricM,3));
        se2plot = zeros(2,size(allMetricM,2),size(allMetricM,3));
        
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
                
            end
        end
        
        if ~exist([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)],'file')
            mkdir([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)])
        end
        save([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}],'lmeTBL','tMap','pMap','y2plot','se2plot','Para');
        
        
    end
    
    
    %% section5: calculate band power %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'BP') | strcmp(calculate,'BPall')
        CondIndexM = ones(size(subIndexM));
        CondIndexS = 2*ones(size(subIndexS));
        allMetricM = squeeze(allMetricM);
        allMetricS = squeeze(allMetricS);
        
        tMap = zeros(1,size(allMetricM,2));
        pMap = zeros(1,size(allMetricM,2));
        y2plot = zeros(2,size(allMetricM,2));
        se2plot = zeros(2,size(allMetricM,2));
        yraw = zeros(2,size(allMetricM,2));
        seraw = zeros(2,size(allMetricM,2));
        
        strlen = 0;
        % TFR point wise LME
        
        for itime = 1:size(allMetricM,2)
            
            s = ['Calculating tf point: ' num2str(itime) '/' num2str(size(allMetricS,2)) 'times'];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;
            
            frameDataM = double(squeeze(allMetricM(:,itime)));
            % skip nan point
            if ~any(frameDataM,'all')
                continue
            end
            frameDataS = double(squeeze(allMetricS(:,itime)));
            
            
            %             lmeTBL = table([frameDataM;frameDataS],[CondIndexM;CondIndexS],[subIndexM;subIndexS], ...
            %                 [elecIndexM;elecIndexS],[trlIndexM;trlIndexS], 'VariableNames',{'Y','Cond','Sub','Elec','Trl'});
            %             lmeTBL.Cond = nominal(lmeTBL.Cond);
            %             lmeTBL.Sub = nominal(lmeTBL.Sub);
            %             lmeTBL.Elec = nominal(lmeTBL.Elec);
            %             lmeTBL.Trl = nominal(lmeTBL.Trl);
            %             lmeStruct = fitlme(lmeTBL,'Y~Cond+(1|Sub)+(1|Elec)','fitmethod','reml','DummyVarCoding','effects');
            
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
            y2plot(:,itime) = [nanmean(obsVal1);nanmean(obsVal2)];
            se2plot(:,itime) = [nanstd(obsVal1)./sqrt(numel(obsVal1));nanstd(obsVal2)./sqrt(numel(obsVal2))];
            rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
            rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
            yraw(:,itime) = [mean(rawVal1);mean(rawVal2)];
            seraw(:,itime) = [std(rawVal1)./sqrt(numel(rawVal1));std(rawVal2)./sqrt(numel(rawVal2))];
        end
        
        if ~exist([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)],'file')
            mkdir([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)])
        end
        save([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}], ...
            'lmeTBL','tMap','pMap','y2plot','se2plot','yraw','seraw','Para');
        
        
    end
    
    
    
    %% section6: calculate PAC %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'PAC') | strcmp(calculate,'PACold')
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
                y2plot(:,ifreq,itime) = [nanmean(obsVal1);nanmean(obsVal2)];
                se2plot(:,ifreq,itime) = [nanstd(obsVal1)./sqrt(numel(obsVal1));nanstd(obsVal2)./sqrt(numel(obsVal2))];
                rawVal1 = lmeTBL.Y(lmeTBL.Cond=='1');
                rawVal2 = lmeTBL.Y(lmeTBL.Cond=='2');
                yraw(:,ifreq,itime) = [nanmean(rawVal1);nanmean(rawVal2)];
                seraw(:,ifreq,itime) = [nanstd(rawVal1)./sqrt(numel(rawVal1));nanstd(rawVal2)./sqrt(numel(rawVal2))];
            end
        end
        
        if ~exist([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)],'file')
            mkdir([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4)])
        end
        save([resultPath calculate filesep ROIAtlas{iatlas}(1:end-4) filesep ROIText{iatlas}], ...
            'lmeTBL','tMap','pMap','y2plot','se2plot','yraw','seraw','Para');
        
    end
    
    
end

varargout{1} = toc/3600;


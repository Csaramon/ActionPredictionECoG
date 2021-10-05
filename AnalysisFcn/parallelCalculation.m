function varargout = parallelCalculation(calculate)

tic;
if nargin < 1
    calculate = 'MPD'
end

% initialize base path and toolbox
if strcmpi(computer,'PCWIN64')
    addpath('C:\Users\qin2\Documents\MATLAB\toolbox')
    addpath('C:\Users\qin2\Documents\MATLAB\toolbox\fieldtrip-20210418')
    basePath = 'C:\Users\qin2\Documents\ActionPrediction\';
    
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

for isub = 1:numel(allsub)
    subname = allsub{isub};
    subPath = [basePath filesep 'Data' filesep subname filesep];
    dataPath = [subPath filesep 'Analysis' filesep];
    fprintf(['\n Currently calculating subject: ' subname])
    
    %% section1: calculate ERP %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'ERP')
        
        statP = 0.05; % significance level for statistical test
        clusterP = 0.05; % significance level for cluster
        %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         if exist([dataPath subname 'CommonAverage_trlData.mat'],'file')
        %             a = load([dataPath subname 'CommonAverage_trlData']);
        %         else
        %             a = load([dataPath subname 'CommonAverageAndBipolar_trlData']);
        %         end
        if exist([dataPath subname 'LAR_preprocData.mat'],'file')
            a = load([dataPath subname 'LAR_preprocData']);
        end
        c = fieldnames(a);
        trlData = a.(c{1});
        
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
        
        
        % datafolder setting
        datafolder = ['ERPLAR' filesep subname filesep];
        if ~exist([resultPath datafolder],'file')
            mkdir(resultPath,datafolder)
        end
        
        cfg = [];
        cfg.channel = 'all';
        cfg.keeptrials = 'yes';
        timelockM = ft_timelockanalysis(cfg,trlDataM);
        timelockS = ft_timelockanalysis(cfg,trlDataS);
        
        maskFDR = zeros(numel(timelockM.label),numel(timelockM.time));
        for ichan = 1:numel(trlData.label)
            tiemlockMdata = squeeze(timelockM.trial(:,ichan,:));
            tiemlockSdata = squeeze(timelockS.trial(:,ichan,:));
            
            [~,p] =  ttest2(tiemlockMdata,tiemlockSdata,'alpha',statP);
            [~,Q] = myfdr(p);
            chanMask = Q;
            chanMask(Q<statP) = 1;
            chanMask(Q>=statP) = 0;
            maskFDR(ichan,:) = chanMask;
            
        end
        
        
        %%%%%%%%-------- Cluster based Permutation Test --------%%%%%%%%
        
        cfg = [];
        cfg.method            = 'montecarlo';           % use the Monte Carlo Method to calculate the significance probability
        cfg.statistic         = 'indepsamplesT';        % use the independent samples T-statistic as a measure to evaluate the effect at the sample level
        cfg.correctm          = 'cluster';
        cfg.clusteralpha      = clusterP;                   % alpha level of the sample-specific test statistic that will be used for thresholding
        cfg.clustertail       = 0;
        cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
        cfg.tail              = 0;                      % -1, 1 or 0 (default = 0); one-sided or two-sided test
        cfg.correcttail       = 'prob';                 % the two-sided test implies that we do non-parametric two tests
        cfg.alpha             = statP;                   % alpha level of the permutation test
        cfg.numrandomization  = 1000;                   % number of draws from the permutation distribution
        cfg.design  = [1*ones(size(timelockM.trial,1),1);2*ones(size(timelockS.trial,1),1)];
        % cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
        cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
        cfg.channel           = 'all';
        cfg.neighbours        = [];
        
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
        
        
        timelockStat = ft_timelockstatistics(cfg, timelockM, timelockS);
        
        save([resultPath datafolder filesep 'TimelockData'],'timelockM','timelockS','timelockStat','maskFDR');
        
        
        %         % plot significance
        %         for ichan = 1:numel(trlData.label)
        %             tiemlockMdata = squeeze(timelockM.trial(:,ichan,:));
        %             tiemlockSdata = squeeze(timelockS.trial(:,ichan,:));
        %             % plot erp
        %             ff = figure('Name', [subname ' ERP in ' timelockM.label{ichan}], 'NumberTitle', 'off');
        %             subplot(2,1,1)
        %             hlineM = shadedErrorBar(timelockM.time, tiemlockMdata, {@mean, @(x) std(x)/sqrt(size(tiemlockMdata,1))},{'color',[255,106,106]/255},1);
        %             hold on
        %             hlineS = shadedErrorBar(timelockS.time, tiemlockSdata, {@mean, @(x) std(x)/sqrt(size(tiemlockSdata,1))},{'color',[30,144,255]/255},1);
        %             ylim([-50 50])
        %             legend([hlineM.mainLine,hlineS.mainLine],['Matched'],['Scrambled'])
        %             [~,p] =  ttest2(tiemlockMdata,tiemlockSdata,'alpha',statP);
        %             [~,Q] = myfdr(p);
        %             subplot(2,1,2)
        %             hP = plot(timelockM.time,p);
        %             hold on
        %             sigh = plot(timelockM.time,repmat(statP,size(timelockM.time)));
        %             sighfdr = plot(timelockM.time,Q);
        %             ylim([0 1])
        %
        %             saveas(ff,[resultPath datafolder filesep timelockM.label{ichan}]);
        %
        %             close(ff)
        %
        %         end
        
    end
    
    
    
    %% section2: calculate ITC %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'ITC')
        
        statP = 0.05; % significance level for statistical test
        clusterP = 0.05; % significance level for cluster
        toi = -0.5:0.05:0.8;
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
        
        
        % Region of Interest include:
        ROIIndex = {[9:12,76:77],[5:6],[15,16,21:22,27:28],[23:24,29:30],[35,36],[1:2],[55],[62:65,100:101,102:103],[13,14],[1:4]}; %62:65,100:101,102:103
        ROIAtlas = {'natlas.nii','natlas.nii','natlas.nii','natlas.nii','natlas.nii','hMTd1.nii','natlas.nii','natlas.nii','natlas.nii','perc_VTPM.nii.gz'};
        ROIText = {'IFG','SFG','aTC','pTC','SPL','hMT','ACC','MTL','preCentral','V1'};
        roiDist = 2; % maximum distance between electrodes and ROI voxels
        
        %         for iatlas = 1:numel(ROIIndex)
        %
        %             cfg = [];
        %             % load altlas infomation
        %             aparc_nii = load_nifti([basePath ROIAtlas{iatlas}]);
        %             aparc_vol = round(aparc_nii.vol);
        %             allroi_index = [];
        %             %%%------ For maxprobabilistic map ------%%%
        %             if ndims(aparc_vol) == 3
        %                 for iroi = ROIIndex{iatlas}
        %                     temp_indx = find(aparc_vol==iroi);
        %                     allroi_index = [allroi_index;temp_indx];
        %                 end
        %             elseif ndims(aparc_vol) == 4
        %                 %%%------ For probabilistic map ------%%%
        %                 aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{iatlas}),4));
        %                 allroi_index = find(aparc_vol>=10);  % minimum probability
        %             end
        
        % atlas parcellation coordinates
        %             [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        %             ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        %             ras_coordiantes = ras_coordiantes(1:3,:)';
        %             aparc_coordiantes = ras_coordiantes;
        
        % load MNI coordinates
        %             elecpos_raw = importdata([subPath filesep 'brain3D/MNI152_coordinates_ras.txt']);
        %             arraynum = tabulate(round(elecpos_raw(:,4)/100));
        %             elecrm = cumsum(arraynum(:,2));
        %             elecpos_raw(elecrm,:) = [];
        %
        % choose ROI electrodes
        %             tempdev = pdist2(elecpos_raw(:,1:3),aparc_coordiantes);
        %             [it,jt] = find(tempdev <=roiDist);
        %             ROIelec = unique(it);
        
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
                    
                    [~,p] = ttest2(CorrM(ichan,it,:) ,CorrS(ichan,it,:),'alpha',statP);
                    
                    CorrP(ichan,it) = p;
                    maskP(ichan,it) = p < (0.05/numel(toi));
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
    
    
    %% section2-2: calculate Inter-view Correlation %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'IVC')
        
        %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist([dataPath subname 'LAR_rerefData.mat'],'file')
            a = load([dataPath subname 'LAR_rerefData']);
        end
        c = fieldnames(a);
        rerefData = a.(c{1});
        
        if exist([dataPath subname '_eventdata.mat'],'file')
            a = load([dataPath subname '_eventdata']);
        end
        c = fieldnames(a);
        camInfo = a.(c{1});
        
        load([dataPath subname  'IVC']); 
        
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
        
        % theta band
        % time-frequency decomposition (wavelet)
%         cfg              = [];
%         cfg.output       = 'pow';
%         cfg.method       = 'wavelet';
%         cfg.foi          = foi(foi>3 & foi<8);
%         cfg.width        =  width(foi>3 & foi<8);
%         cfg.toi          = 'all';
%         cfg.precision = 'single';
%         
%         ft_warning off
%         freq = ft_freqanalysis(cfg,rerefData);
%         
%         allTS = freq.powspctrm;
%         clear freq
%         
%         % calculate IVC
%         InterViewCorr
%         IVC.Intact.theta = IVC_I;
%         IVC.Intact.thetaShf = IVC_Ishf;
%         IVC.Scamble.theta = IVC_S;
%         IVC.Scamble.thetaShf = IVC_Sshf;
%         IVC.Intact_Scamble.theta = IVC_I_S;
%         IVC.Intact_Scamble.thetaShf = IVC_I_Sshf;
%         IVC.Scamble_Intact.theta = IVC_S_I;
%         IVC.Scamble_Intact.thetaShf = IVC_S_Ishf;
        
        % alpha band
        % time-frequency decomposition (wavelet)
%         cfg              = [];
%         cfg.output       = 'pow';
%         cfg.method       = 'wavelet';
%         cfg.foi          = foi(foi>8 & foi<13);
%         cfg.width        =  width(foi>8 & foi<13);
%         cfg.toi          = 'all';
%         cfg.precision = 'single';
%         
%         ft_warning off
%         freq = ft_freqanalysis(cfg,rerefData);
%         
%         allTS = freq.powspctrm;
%         clear freq
%         
%         % calculate IVC
%         InterViewCorr
%         IVC.Intact.alpha = IVC_I;
%         IVC.Intact.alphaShf = IVC_Ishf;
%         IVC.Scamble.alpha = IVC_S;
%         IVC.Scamble.alphaShf = IVC_Sshf;
%         IVC.Intact_Scamble.alpha = IVC_I_S;
%         IVC.Intact_Scamble.alphaShf = IVC_I_Sshf;
%         IVC.Scamble_Intact.alpha = IVC_S_I;
%         IVC.Scamble_Intact.alphaShf = IVC_S_Ishf;
        
          % low beta band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>13 & foi<20);
        cfg.width        =  width(foi>13 & foi<20);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);
        
        allTS = freq.powspctrm;
        clear freq
        
        % calculate IVC
        InterViewCorr
        IVC.Intact.lbeta = IVC_I;
        IVC.Intact.lbetaShf = IVC_Ishf;
        IVC.Scamble.lbeta = IVC_S;
        IVC.Scamble.lbetaShf = IVC_Sshf;
        IVC.Intact_Scamble.lbeta = IVC_I_S;
        IVC.Intact_Scamble.lbetaShf = IVC_I_Sshf;
        IVC.Scamble_Intact.lbeta = IVC_S_I;
        IVC.Scamble_Intact.lbetaShf = IVC_S_Ishf;        

        % high beta band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>20 & foi<30);
        cfg.width        =  width(foi>20 & foi<30);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);
        
        allTS = freq.powspctrm;
        clear freq
        
        % calculate IVC
        InterViewCorr
        IVC.Intact.hbeta = IVC_I;
        IVC.Intact.hbetaShf = IVC_Ishf;
        IVC.Scamble.hbeta = IVC_S;
        IVC.Scamble.hbetaShf = IVC_Sshf;
        IVC.Intact_Scamble.hbeta = IVC_I_S;
        IVC.Intact_Scamble.hbetaShf = IVC_I_Sshf;
        IVC.Scamble_Intact.hbeta = IVC_S_I;
        IVC.Scamble_Intact.hbetaShf = IVC_S_Ishf;


        % beta band
        % time-frequency decomposition (wavelet)
%         cfg              = [];
%         cfg.output       = 'pow';
%         cfg.method       = 'wavelet';
%         cfg.foi          = foi(foi>13 & foi<30);
%         cfg.width        =  width(foi>13 & foi<30);
%         cfg.toi          = 'all';
%         cfg.precision = 'single';
%         
%         ft_warning off
%         freq = ft_freqanalysis(cfg,rerefData);
%         
%         allTS = freq.powspctrm;
%         clear freq
%         
%         % calculate IVC
%         InterViewCorr
%         IVC.Intact.beta = IVC_I;
%         IVC.Intact.betaShf = IVC_Ishf;
%         IVC.Scamble.beta = IVC_S;
%         IVC.Scamble.betaShf = IVC_Sshf;
%         IVC.Intact_Scamble.beta = IVC_I_S;
%         IVC.Intact_Scamble.betaShf = IVC_I_Sshf;
%         IVC.Scamble_Intact.beta = IVC_S_I;
%         IVC.Scamble_Intact.betaShf = IVC_S_Ishf;
        
        
        
        % lowgamma band
        % time-frequency decomposition (wavelet)
%         cfg              = [];
%         cfg.output       = 'pow';
%         cfg.method       = 'wavelet';
%         cfg.foi          = foi(foi>30 & foi<60);
%         cfg.width        =  width(foi>30 & foi<60);
%         cfg.toi          = 'all';
%         cfg.precision = 'single';
%         
%         ft_warning off
%         freq = ft_freqanalysis(cfg,rerefData);
%         
%         allTS = freq.powspctrm;
%         clear freq
%         
%         % calculate IVC
%         InterViewCorr
%         IVC.Intact.lgamma = IVC_I;
%         IVC.Intact.lgammaShf = IVC_Ishf;
%         IVC.Scamble.lgamma = IVC_S;
%         IVC.Scamble.lgammaShf = IVC_Sshf;
%         IVC.Intact_Scamble.lgamma = IVC_I_S;
%         IVC.Intact_Scamble.lgammaShf = IVC_I_Sshf;
%         IVC.Scamble_Intact.lgamma = IVC_S_I;
%         IVC.Scamble_Intact.lgammaShf = IVC_S_Ishf;
        
        
        % highgamma band
        % time-frequency decomposition (wavelet)
%         cfg              = [];
%         cfg.output       = 'pow';
%         cfg.method       = 'wavelet';
%         cfg.foi          = foi(foi>60 & foi<120);
%         cfg.width        =  width(foi>60 & foi<120);
%         cfg.toi          = 'all';
%         cfg.precision = 'single';
%         
%         ft_warning off
%         freq = ft_freqanalysis(cfg,rerefData);
%         
%         allTS = freq.powspctrm;
%         clear freq
%         
%         % calculate IVC
%         InterViewCorr
%         IVC.Intact.hgamma = IVC_I;
%         IVC.Intact.hgammaShf = IVC_Ishf;
%         IVC.Scamble.hgamma = IVC_S;
%         IVC.Scamble.hgammaShf = IVC_Sshf;
%         IVC.Intact_Scamble.hgamma = IVC_I_S;
%         IVC.Intact_Scamble.hgammaShf = IVC_I_Sshf;
%         IVC.Scamble_Intact.hgamma = IVC_S_I;
%         IVC.Scamble_Intact.hgammaShf = IVC_S_Ishf;
        
        
        % save data
        save([dataPath subname  'IVC'],'IVC');
        
        
    end
    
    
    
    %% section3: calculate power spectrum and peak frequency %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'PeakFrequency')
        % set parameter
        freqRange = [1 128];
        timeRange = [0 1];
        timeTFR = 0.02;
        
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
            cfg.toi          = min(trlData.time{1}):timeTFR:max(trlData.time{1}); %'all';
            
            ft_warning off
            freq         = ft_freqanalysis(cfg,rerefSingle);
            
            % select data
            toi = cfg.toi>=timeRange(1) & cfg.toi<=timeRange(2);
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
    
    
    %% section3.2: plot powerspectrum and peak in MNI coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'PlotPowerPeak')
        directorynameall{1} = [resultPath 'PowerSpectrum/1-128Hz'];
        % directorynameall{2} = '/Users/apple/Desktop/WorkFLow/20190427PowerSpectrum/replay_TrialSeperate';
        % directorynameall{3} = '/Users/apple/Desktop/WorkFLow/20190427PowerSpectrum/bistable_WholeTime';
        for idd = 1%1:numel(directorynameall)
            
            
            directoryname = directorynameall{idd};
            D=dir([directoryname filesep '*.fig']);
            listname = char({D.name});
            
            
            % plot MNI152 surface
            vox2ras = [-1,0,0,127;0,0,1,-145;0,-1,0,147;0,0,0,1];
            ras_tkr2vox = [-1,0,0,128;0,0,-1,128;0,1,0,128;0,0,0,1];
            [~,MNI152path] = unix('echo $SUBJECTS_DIR ');
            MNI152path = [MNI152path(1:end-1) '/cvs_avg35_inMNI152/'];
            
            % Region of Interest include:
            ROIIndex = {[61,62],[63,64],[1]};
            ROIAtlas = {'fsAnatomyMacro.nii','fsAnatomyMacro.nii','T1mask.nii'};
            ROIText = {'InferiorParietalLobe','SupraMarginal','WholeBrain'};
            roiDist = 1; % maximum distance between electrodes and ROI voxels
            
            peak_frequency_roi = [];
            for iatlas = 1:numel(ROIAtlas)
                
                aparc_nii = load_nifti([basePath '/Atlas/' ROIAtlas{iatlas}]);
                aparc_vol = round(aparc_nii.vol);
                
                allroi_index = [];
                if ndims(aparc_vol) == 3
                    %%%------ For maxprobabilistic map ------%%%
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
                
                % set figure handle
                stereoview = figure('Name','Distribution of Peak Frequencys','NumberTitle','of');
                stereoaxes = axes('Parent',stereoview);
                hold(stereoaxes,'on')
                axis(stereoaxes,'equal');
                axis(stereoaxes,'off');
                % plot pial surface
                [vertex_coords1, faces1]=read_surf([MNI152path '/surf/lh.pial']);
                faces1 = faces1 + 1;
                [vertex_coords2, faces2]=read_surf([MNI152path '/surf/rh.pial']);
                faces2 = faces2 + 1;
                stereonum = [size(vertex_coords1,1) size(faces1,1)];
                faces = [faces1;faces2+size(vertex_coords1,1)];
                vertices = [vertex_coords1;vertex_coords2];
                
                vertices = vox2ras*ras_tkr2vox*[vertices ones(size(vertices,1),1)]';
                vertices = vertices(1:3,:)';
                stereo = patch(struct(...
                    'vertices', vertices, 'faces', faces), ...
                    'Parent',gca, ...
                    'FaceColor',[230,228,216]./255, ...
                    'FaceAlpha',0.25, ...
                    'EdgeColor', 'none');
                
                light = camlight;
                set(light,'position',[0 0 2000])
                set(light,'position',get(gca,'cameraPosition'))
                lighting gouraud;
                material dull
                
                
                filefolder = [directoryname filesep ROIText{iatlas}];
                if ~exist(filefolder,'file')
                    mkdir(filefolder)
                end
                
                mycolormap = jet(64);
                peak_frequency_all = [];
                activecoordinates_all = [];
                numElecFoi = 0;
                numElecAll = 0;
                numElecOsc = 0;
                
                
                fpow = figure('Name',[ROIText{iatlas} '_powerspctrum']);
                hfp = axes('Parent',fpow);
                hold on
                
                
                for ielec = 1:size(listname,1)
                    tempHF=  openfig([directoryname filesep deblank(listname(ielec,:))],'invisible');
                    
                    activecoordinates =      getappdata(tempHF,'MNICoordinate');
                    spec_to_plot =  getappdata(tempHF,'power');
                    foi =    getappdata(tempHF,'frequency');
                    peak_frequency =   getappdata(tempHF,'peakFrequency');
                    
                    
                    hp = plot(hfp,log2(foi),log2(spec_to_plot)./sum(log2(spec_to_plot)));
                    set(hfp,'xtick',[0:7],'xticklabel',[1,2,4,8,16,32,64,128]);
                    xlim([0 7])
                    
                    tempdev = pdist2(activecoordinates(:,1:3),aparc_coordiantes);
                    if ~any(tempdev <=roiDist) % exclude coordinates out of range
                        close(tempHF)
                        continue
                    else
                        numElecAll = numElecAll+1;
                    end
                    
                    if ~isempty(peak_frequency)
                        numElecOsc = numElecOsc+1;
                        %                 if any(peak_frequency>=2 & peak_frequency<=7)
                        %                     numElecFoi = numElecFoi+1;
                        %                     setappdata(tempHF,'OSCstate',1);
                        %                     saveas(tempHF,[directoryname filesep deblank(listname(ielec,:))]);
                        %                 else
                        %                     setappdata(tempHF,'OSCstate',0);
                        %                     saveas(tempHF,[directoryname filesep deblank(listname(ielec,:))]);
                        %                 end
                        relativeCoordinates = [];
                        for ip = 1:numel(peak_frequency)
                            relativeCoordinates(ip,1:3) = activecoordinates(1,1:3)+min(1,ip-1)*rand(1,3);
                            relativeCoordinates(ip,4) = activecoordinates(1,4);
                        end
                        activecoordinates_all = [activecoordinates_all;relativeCoordinates];
                        peak_frequency_all = [peak_frequency_all,peak_frequency];
                        
                    else
                        setappdata(tempHF,'OSCstate',0);
                        saveas(tempHF,[directoryname filesep deblank(listname(ielec,:))]);
                        
                    end
                    
                    close(tempHF)
                    tempHF = [];
                    
                end
                
                for ie = 1:numel(peak_frequency_all)
                    ecolor = mycolormap(round(peak_frequency_all(ie)/max(peak_frequency_all)*64),:);
                    he(ie) = scatter3(stereoaxes,activecoordinates_all(ie,1),activecoordinates_all(ie,2),activecoordinates_all(ie,3),60,ecolor,'fill');
                    drawnow
                end
                
                if isempty(peak_frequency_all)
                    continue
                end
                
                hh = figure('Name', 'Histogram of peak frequencys','NumberTitle', 'off');
                hha = axes('Parent',hh);
                %     hhist = histogram(hha,peak_frequency_all,'BinMethod','fd');
                hhist = histogram(hha,peak_frequency_all,'BinWidth',1);
                xlabel(hha,'Peak Frequencys(Hz)')
                ylabel(hha,'Counts')
                title(stereoaxes,['Electrode Peak Frequency Distribution (' num2str(numElecFoi) 'in ' num2str(numElecOsc) ' in ' num2str(numElecAll) ')'])
                title(hha,['Peak Frequency Distribution in ' ROIText{iatlas}])
                
                colormap(stereoaxes,mycolormap)
                hcolorbar = colorbar('peer',stereoaxes);
                set(hcolorbar,'position',[0.925 0.25 0.0333 0.55])
                caxis(stereoaxes,[min(peak_frequency_all) max(peak_frequency_all)])
                
                
                
                saveas(hh,[filefolder filesep 'FrequencyPeakHistogram']);
                saveas(stereoview,[filefolder filesep 'FrequencyPeakDistribution']);
                saveas(fpow,[filefolder filesep 'PowerSpectrum']);
                
                close all
                
                
                peak_frequency_roi{1,iatlas} = [peak_frequency_all'];
                
            end
            
            
            hv = violin(peak_frequency_roi,'x',[1:numel(peak_frequency_roi)]);
            set(gca,'XTickLabel',ROIText(1:iatlas));
            save([directoryname filesep 'FreqvuencyPeak'],'peak_frequency_roi');
            saveas(hv,[directoryname filesep 'FreqvuencyPeakViolin']);
            close all
            clearvars -EXCEPT idd directorynameall
        end
    end
    
    
    %% section4: calculate time frequency representation %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'TFR')
        
        timeWinTFR = [-0.5 1]; % time window for time frequency results
        timeTFR = 0.01; % time step for time frequency results
        
        
        %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist([dataPath subname 'LAR_trlData.mat'],'file')
            a = load([dataPath subname 'LAR_trlData']);
        end
        c = fieldnames(a);
        trlData = a.(c{1});
        
        cfg = [];
        cfg.demean = 'yes';
        cfg.detrend = 'yes';
        
        trlData = ft_preprocessing(cfg,trlData);
        
        % seperate conditions
        cfg = [];
        %         cfg.channel = ERPchan{isub};
        cfg.trials = find(trlData.trialinfo(:,1)==0);
        trlDataM = ft_selectdata(cfg,trlData);
        
        cfg = [];
        %         cfg.channel = ERPchan{isub};
        cfg.trials = find(trlData.trialinfo(:,1)==1);
        trlDataS = ft_selectdata(cfg,trlData);
        
        
        % datafolder setting
        datafolder = ['TFR' filesep subname filesep];
        if ~exist([resultPath datafolder],'file')
            mkdir(resultPath,datafolder)
        end
        
        %         wavelet
        cfg              = [];
        cfg.output       = 'fourier';
        cfg.method       = 'wavelet';
        cfg.foi          = logspace(log10(2),log10(128),32);
        cfg.width        =  logspace(log10(3),log10(10),32); % adjustive cycles
        cfg.toi          = min(timeWinTFR):timeTFR:max(timeWinTFR); %'all';
        cfg.keeptrials   = 'yes';
        cfg.precision = 'single';
        cfg.pad = 'nextpow2';
        
        % STFFT
        %                 cfg              = [];
        %         cfg.output       = 'fourier';
        %         cfg.method       = 'mtmconvol';
        %         cfg.taper        = 'hanning';
        %         cfg.foi          = logspace(log10(2),log10(128),32);
        %         cycles        =  logspace(log10(3),log10(10),32); % adjustive cycles
        %         cfg.t_ftimwin  = cycles./cfg.foi;
        %         cfg.toi          = min(timeWinTFR):timeTFR:max(timeWinTFR); %'all';
        %         cfg.keeptrials   = 'yes';
        %         cfg.precision = 'single';
        
        %%%---------- hilbert ----------%%%
        %         cfg              = [];
        %         cfg.output       = 'fourier';
        %         cfg.method       = 'hilbert';
        %         cfg.foi          = logspace(log10(2),log10(128),32);
        %         cfg.width        =  0.3*cfg.foi;
        %         cfg.toi          = min(timeWinTFR):timeTFR:max(timeWinTFR); %'all';
        %         cfg.keeptrials   = 'yes';
        %         cfg.filttype = 'firws';
        %         cfg.filtorder = nan;
        %         cfg.filtdir = 'onepass-zerophase';
        %         cfg.precision = 'single';
        
        % multi-taper window
        %         cfg              = [];
        %         cfg.output       = 'fourier';
        %         cfg.method       = 'mtmconvol';
        %         cfg.foi          = logspace(log10(2),log10(128),32);
        %         cfg.t_ftimwin  = 5./cfg.foi;
        %         cfg.tapsmofrq  = 0.4 *cfg.foi; % tapers=2*tw*fw-1
        %         cfg.toi          = min(trlDataM.time{1}):timeTFR:max(trlDataM.time{1}); %'all';
        %         cfg.keeptrials   = 'yes';
        %         cfg.precision = 'single';
        
        ft_warning off
        freqM = ft_freqanalysis(cfg,trlDataM);
        freqS = ft_freqanalysis(cfg,trlDataS);
        
        %%%%%%%%-------- Cluster based Permutation Test --------%%%%%%%%
        
        %         cfg = [];
        %         cfg.channel           = 'all';
        %         cfg.method            = 'montecarlo';           % use the Monte Carlo Method to calculate the significance probability(Ref: ft_statistics_montecarlo)
        %         cfg.statistic         = 'indepsamplesT';        % use the independent samples T-statistic as a measure to evaluate the effect at the sample level
        %         cfg.correctm          = 'cluster';
        %         cfg.clusteralpha      = clusterP;                   % alpha level of the sample-specific test statistic that will be used for thresholding
        %         cfg.clustertail       = 0;
        %         cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
        %         cfg.clusterthreshold = 'nonparametric_individual';
        %         cfg.tail              = 0;                      % -1, 1 or 0 (default = 0); one-sided or two-sided test
        %         cfg.alpha             = statP;                   % alpha level of the permutation test
        %         cfg.numrandomization  = 1000;                   % number of draws from the permutation distribution
        %
        %         design = zeros(1,size(freqM.powspctrm,1) + size(freqS.powspctrm,1));
        %         design(1,1:size(freqM.powspctrm,1)) = 1;
        %         design(1,(size(freqM.powspctrm,1)+1):(size(freqM.powspctrm,1)+...
        %             size(freqS.powspctrm,1))) = 2;
        %         cfg.design            = design; % design matrix, note the transpose
        %         cfg.ivar              = 1;                      % the index of the independent variable in the design matrix
        %         cfg.neighbours        = [];                     % there are no spatial neighbours, only in time and frequency
        %         cfg.parameter   = 'powspctrm';
        %
        %         freqStat = ft_freqstatistics(cfg, freqM, freqS);
        %%%%%%%%-------- Cluster based Permutation Test --------%%%%%%%%
        
        
        %         save([resultPath datafolder filesep 'FreqData'],'freqM','freqS','freqStat');
        save([resultPath datafolder filesep 'FreqData'],'freqM','freqS');
        
        
    end
    
    
    
   %% section5: calculate grand power difference %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(calculate,'MPD')
        
        %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist([dataPath subname 'LAR_rerefData.mat'],'file')
            a = load([dataPath subname 'LAR_rerefData']);
        end
        c = fieldnames(a);
        rerefData = a.(c{1});
        
        load([dataPath subname  'PowDiff']);
        
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
        
        
        % theta band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>=3 & foi<=8);
        cfg.width        =  width(foi>=3 & foi<=8);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);
 
                % no filtering
%         allTS = freq.powspctrm;
%         % calculate MPD
%         MeanPowDiff
%         PowDiff.theta = PD_I_S;
%         PowDiff.thetaShf = PD_shf;        
        allTS = freq.powspctrm;
        % low pass filtering
        for ifreq = 1:size(allTS,2)
            [allTS(:,ifreq,:), ~, ~] = ft_preproc_lowpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
        end
        % calculate MPD
        MeanPowDiff
        lpPowDiff.theta = PD_I_S;
        lpPowDiff.thetaShf = PD_shf;
        
        allTS = freq.powspctrm;
        % high pass filtering
        for ifreq = 1:size(allTS,2)
            [allTS(:,ifreq,:), ~, ~] = ft_preproc_highpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
        end
        % calculate MPD
        MeanPowDiff
        hpPowDiff.theta = PD_I_S;
        hpPowDiff.thetaShf = PD_shf;
        
                 % alpha band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>=8 & foi<=13);
        cfg.width        =  width(foi>=8 & foi<=13);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);

                % no filtering
        allTS = freq.powspctrm;
        % calculate MPD
        MeanPowDiff
        PowDiff.alpha = PD_I_S;
        PowDiff.alphaShf = PD_shf;
%         allTS = freq.powspctrm;
%         % low pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_lowpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         lpPowDiff.alpha = PD_I_S;
%         lpPowDiff.alphaShf = PD_shf;
%         
%         allTS = freq.powspctrm;
%         % high pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_highpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         hpPowDiff.alpha = PD_I_S;
%         hpPowDiff.alphaShf = PD_shf;
        
                % low beta band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>=13 & foi<=20);
        cfg.width        =  width(foi>=13 & foi<=20);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);

                % no filtering
        allTS = freq.powspctrm;
        % calculate MPD
        MeanPowDiff
        PowDiff.lbeta = PD_I_S;
        PowDiff.lbetaShf = PD_shf;
%         allTS = freq.powspctrm;
%         % low pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_lowpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         lpPowDiff.lbeta = PD_I_S;
%         lpPowDiff.lbetaShf = PD_shf;
%         
%         allTS = freq.powspctrm;
%         % high pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_highpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         hpPowDiff.lbeta = PD_I_S;
%         hpPowDiff.lbetaShf = PD_shf;


        % high beta band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>=20 & foi<=30);
        cfg.width        =  width(foi>=20 & foi<=30);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);
        
                % no filtering
        allTS = freq.powspctrm;
        % calculate MPD
        MeanPowDiff
        PowDiff.hbeta = PD_I_S;
        PowDiff.hbetaShf = PD_shf;
        
%             allTS = freq.powspctrm;
%         % low pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_lowpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         lpPowDiff.hbeta = PD_I_S;
%         lpPowDiff.hbetaShf = PD_shf;
%         
%         allTS = freq.powspctrm;
%         % high pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_highpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         hpPowDiff.hbeta = PD_I_S;
%         hpPowDiff.hbetaShf = PD_shf;
 
        
                % low gamma band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>=30 & foi<=60);
        cfg.width        =  width(foi>=30 & foi<=60);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);
        
                % no filtering
        allTS = freq.powspctrm;
        % calculate MPD
        MeanPowDiff
        PowDiff.lgamma = PD_I_S;
        PowDiff.lgammaShf = PD_shf;
        
%               allTS = freq.powspctrm;
%         % low pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_lowpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         lpPowDiff.lgamma = PD_I_S;
%         lpPowDiff.lgammaShf = PD_shf;
%         
%         allTS = freq.powspctrm;
%         % high pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_highpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         hpPowDiff.lgamma = PD_I_S;
%         hpPowDiff.lgammaShf = PD_shf;
        
        
        % high gamma band
        % time-frequency decomposition (wavelet)
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'wavelet';
        cfg.foi          = foi(foi>=60 & foi<=120);
        cfg.width        =  width(foi>=60 & foi<=120);
        cfg.toi          = 'all';
        cfg.precision = 'single';
        
        ft_warning off
        freq = ft_freqanalysis(cfg,rerefData);
        
        % no filtering
        allTS = freq.powspctrm;
        % calculate MPD
        MeanPowDiff
        PowDiff.hgamma = PD_I_S;
        PowDiff.hgammaShf = PD_shf;
        
%               allTS = freq.powspctrm;
             %         % low pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_lowpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         lpPowDiff.hgamma = PD_I_S;
%         lpPowDiff.hgammaShf = PD_shf;
%         
%         allTS = freq.powspctrm;
%         % high pass filtering
%         for ifreq = 1:size(allTS,2)
%             [allTS(:,ifreq,:), ~, ~] = ft_preproc_highpassfilter(squeeze(allTS(:,ifreq,:)),500,0.5,[],'firws');
%         end
%         % calculate MPD
%         MeanPowDiff
%         hpPowDiff.hgamma = PD_I_S;
%         hpPowDiff.hgammaShf = PD_shf;
        
        % save data
        save([dataPath subname  'PowDiff'],'PowDiff','lpPowDiff','hpPowDiff');
        
        
    end
    
    
end

varargout{1} = toc/3600;
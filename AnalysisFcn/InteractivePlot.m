%%%%%%%%-------- Script for Interactive Plotting --------%%%%%%%%

ft_defaults
subname = 'Patient13';
calculate = 'ERP';


basePath = 'C:\Users\qin2\Documents\ActionPredictionECoG\';
subPath = [basePath subname filesep];
resultPath = [basePath 'Results/'];

if exist([basePath 'Data\' subname '\Analysis\' subname 'LARER_trlData.mat'],'file')
    a = load([basePath 'Data\' subname '\Analysis\' subname 'LARER_trlData.mat']);
end
c = fieldnames(a);
trlData = a.(c{1});
            
%% preprocessing
cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 100;
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
            timelockM = ft_timelockanalysis(cfg,trlDataM);
            timelockS = ft_timelockanalysis(cfg,trlDataS);
%% prepare headshape and layout
pial_lh = ft_read_headshape({[basePath '\Data\' subname '/FsRecon/surf/lh.pial'], ...
    [basePath '\Data\' subname '/FsRecon/surf/rh.pial']});
pial_lh.coordsys = 'acpc';
cfg = [];
cfg.headshape = pial_lh;
cfg.projection = 'orthographic';
% cfg.channel = {'A*'};
cfg.viewpoint = 'left'; % left, right, superior, inferior
% cfg.viewpoint = 'inferior'; % left, right, superior, inferior
cfg.mask = 'convex';
cfg.boxchannel = {'A2', 'A4'};

 %%%%-------- for ERP results --------%%%%
if strncmp(calculate,'ERP',3)
    
    lay = ft_prepare_layout(cfg,timelockM);
    % average the data and plot
%     cfg = [];
%     cfg.channel = 'all';
%     timelockM = ft_timelockanalysis(cfg,timelockM);
%     timelockS = ft_timelockanalysis(cfg,timelockS);
    
    % interactive plotting
    cfg = [];
    cfg.layout = lay;
    cfg.showoutline = 'yes';
    cfg.linewidth = 1.5;
%     cfg.maskparameter = 'mask';
%     cfg.maskfacealpha = 0.5;
%     cfg.channel = 1:104;
    % plot with cluster permutation
%     timelockM.mask = timelockStat.mask;
    
    ft_multiplotER(cfg,timelockM, timelockS)
    
    % plot with FDR correction
    % timelockM.mask = maskFDR;
    %
    % ft_multiplotER(cfg,timelockM, timelockS)
    
        
    %%%%-------- for TFR results --------%%%%
elseif strncmp(calculate,'TFR',3)
    
    lay = ft_prepare_layout(cfg,freqM);
    
    % average the data 
    cfg = [];
    cfg.avgoverrpt = 'yes';
    freqAvgM= ft_selectdata(cfg,freqM);
    freqAvgS= ft_selectdata(cfg,freqS);
    
    cfg = [];
    cfg.parameter    = 'powspctrm';
    cfg.operation    = '(x1-x2) / (x1+x2)';
    
    freqDiff = ft_math(cfg, freqAvgM, freqAvgS);
    freqDiff.mask = freqStat.mask;
    
    cfg = [];
    cfg.layout = lay;
    cfg.showoutline = 'yes';
    cfg.linewidth = 1.5;
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
    
    ft_multiplotTFR(cfg,freqDiff)
    
    
end
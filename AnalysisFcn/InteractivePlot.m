%%%%%%%%-------- Script for Interactive Plotting --------%%%%%%%%

ft_defaults
subname = 'Patient8';
calculate = 'ERP';


basePath = '/Users/qinchaoyi/Desktop/ActionPrediction/';
subPath = [basePath subname filesep];
resultPath = [basePath 'Results/'];

datafile = dir([resultPath calculate filesep subname filesep '*.mat']);
load([datafile.folder filesep datafile.name]);
pial_lh = ft_read_headshape({[basePath subname '/FsRecon/surf/lh.pial'], ...
    [basePath subname '/FsRecon/surf/rh.pial']});
pial_lh.coordsys = 'acpc';


%% prepare headshape and layout
cfg = [];
cfg.headshape = pial_lh;
cfg.projection = 'orthographic';
cfg.channel = {'A*'};
cfg.viewpoint = 'left'; % left, right, superior, inferior
% cfg.viewpoint = 'inferior'; % left, right, superior, inferior
cfg.mask = 'convex';
cfg.boxchannel = {'A1', 'A3'};

 %%%%-------- for ERP results --------%%%%
if strncmp(calculate,'ERP',3)
    
    lay = ft_prepare_layout(cfg,timelockM);
    % average the data and plot
    cfg = [];
    cfg.channel = 'all';
    timelockM = ft_timelockanalysis(cfg,timelockM);
    timelockS = ft_timelockanalysis(cfg,timelockS);
    
    % interactive plotting
    cfg = [];
    cfg.layout = lay;
    cfg.showoutline = 'yes';
    cfg.linewidth = 1.5;
    cfg.maskparameter = 'mask';
    cfg.maskfacealpha = 0.5;
%     cfg.channel = 1:104;
    % plot with cluster permutation
    timelockM.mask = timelockStat.mask;
    
    ft_multiplotER(cfg,timelockM, timelockS)
    
    % plot with FDR correction
    % timelockM.mask = maskFDR;
    %
    % ft_multiplotER(cfg,timelockM, timelockS)
    
        % plot ITC

        datafile = dir([resultPath 'ITC' filesep subname filesep '*.mat']);
        load([datafile.folder filesep datafile.name]);
    cfg = [];
    cfg.layout = lay;
    cfg.showoutline = 'yes';
    cfg.linewidth = 1.5;
%     cfg.channel = 1:104;
        ft_multiplotER(cfg,ITCM, ITCS)
        
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
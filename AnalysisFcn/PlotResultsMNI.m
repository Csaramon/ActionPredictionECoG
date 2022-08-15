%% data path definition
basePath = '/Users/qinchaoyi/Documents/ActionPredictionECoG/';
dataPath = [basePath 'Data/'];
resultPath = [basePath 'Results/'];

%% plot IVC electrodes on MNI brain

allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub

for isub = [1,3,4,6,7,8,9,10]%1:numel(allsub)%[1,3,4,6,7,8,9,10]%1:numel(allsub)
    subname = allsub{isub};
    subPath = [basePath subname filesep];
    dataPath = [subPath filesep 'Analysis' filesep];
    fprintf(['\n Currently calculating subject: ' subname])
    
    
    p=0.05; % threshold for IVC
    if exist([dataPath subname 'IVC.mat'],'file')
        load([dataPath subname 'IVC']);
    end
    
    respElecI = find(IVC.Intact.theta(:,2)<p);
    respElecS = find(IVC.Scamble.theta(:,2)<p);
    respElecInS = find(IVC.Intact.theta(:,2)<p & IVC.Scamble.theta(:,2)<p);
    %         respElecI = find(IVC.Intact.alpha(:,2)<p);
    %     respElecS = find(IVC.Scamble.alpha(:,2)<p);
    %     respElecInS = find(IVC.Intact.alpha(:,2)<p & IVC.Scamble.alpha(:,2)<p);
    %             respElecI = find(IVC.Intact.beta(:,2)<p);
    %     respElecS = find(IVC.Scamble.beta(:,2)<p);
    %     respElecInS = find(IVC.Intact.beta(:,2)<p & IVC.Scamble.beta(:,2)<p);
    %                 respElecI = find(IVC.Intact.lgamma(:,2)<p);
    %     respElecS = find(IVC.Scamble.lgamma(:,2)<p);
    %     respElecInS = find(IVC.Intact.lgamma(:,2)<p & IVC.Scamble.lgamma(:,2)<p);
    %                     respElecI = find(IVC.Intact.hgamma(:,2)<p);
    %     respElecS = find(IVC.Scamble.hgamma(:,2)<p);
    %     respElecInS = find(IVC.Intact.hgamma(:,2)<p & IVC.Scamble.hgamma(:,2)<p);
    
    % load electrodes coordinates
    load([subPath 'FsRecon/brain3D/ft_elec.mat'])
    % project electrode to the left hemisphere
    elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1) = elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1).*-1;
    
    for ie = respElecI'
        
        hI(ie) = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[255 0 0]./255,'fill');
    end
    
    for ie = respElecS'
        
        hS(ie) = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[0 0 255]./255,'fill');
    end
    
    for ie = respElecInS'
        
        hInS(ie) = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[0 255 0]./255,'fill');
    end
    
end

legend([hI(respElecI(1)),hS(respElecS(1)),hInS(respElecInS(1))],['Intact > 0'],['Scrambled > 0'],['Intact & Scrambled > 0'])


%% plot IVC preference electrodes on MNI brain

PlotMNIBrain
allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub

for isub = [1,3,4,6,7,8,9,10]%1:numel(allsub)%[1,3,4,6,7,8,9,10]%1:numel(allsub)
    subname = allsub{isub};
    subPath = [basePath subname filesep];
    dataPath = [subPath filesep 'Analysis' filesep];
    fprintf(['\n Currently calculating subject: ' subname])
    
    
    p=0.05; % threshold for IVC
    if exist([dataPath subname 'IVC.mat'],'file')
        load([dataPath subname 'IVC']);
    end
    %     load([resultPath '/avgpow/fsAnatomyMacro20_30Hz' filesep subname])
    %
    
    %     respElecInS = find(IVC.Intact.theta(:,2)<p & IVC.Intact_Scamble.theta(:,2)<p);
    %         respElecSnI = find(IVC.Scamble.theta(:,2)<p & IVC.Scamble_Intact.theta(:,2)<p);
    % respElecN = setdiff(1:size(IVC.Intact.theta,1),[respElecInS;respElecSnI]);
    
    %     respElecInS = find(IVC.Intact.alpha(:,2)<p & IVC.Intact_Scamble.alpha(:,2)<p);
    %         respElecSnI = find(IVC.Scamble.alpha(:,2)<p & IVC.Scamble_Intact.alpha(:,2)<p);
    % respElecN = setdiff(1:size(IVC.Intact.alpha,1),[respElecInS;respElecSnI]);
    
    %             respElecInS = find(IVC.Intact.lbeta(:,2)<p & IVC.Intact_Scamble.lbeta(:,2)<p);
    %         respElecSnI = find(IVC.Scamble.lbeta(:,2)<p & IVC.Scamble_Intact.lbeta(:,2)<p);
    %         respElecN = setdiff(1:size(IVC.Intact.lbeta,1),[respElecInS;respElecSnI]);
    %             respElecInS = find(IVC.Intact.hbeta(:,2)<p & IVC.Intact_Scamble.hbeta(:,2)<p);
    %         respElecSnI = find(IVC.Scamble.hbeta(:,2)<p & IVC.Scamble_Intact.hbeta(:,2)<p);
    %         respElecN = setdiff(1:size(IVC.Intact.hbeta,1),[respElecInS;respElecSnI]);
    respElecInS = find(IVC.Intact.beta(:,2)<p & IVC.Intact_Scamble.beta(:,2)<p);
    respElecSnI = find(IVC.Scamble.beta(:,2)<p & IVC.Scamble_Intact.beta(:,2)<p);
    respElecN = setdiff(1:size(IVC.Intact.beta,1),[respElecInS;respElecSnI]);
    
    %            respElecInS = find(IVC.Intact.lgamma(:,2)<p & IVC.Intact_Scamble.lgamma(:,2)<p);
    %         respElecSnI = find(IVC.Scamble.lgamma(:,2)<p & IVC.Scamble_Intact.lgamma(:,2)<p);
    %         respElecN = setdiff(1:size(IVC.Intact.lgamma,1),[respElecInS;respElecSnI]);
    %            respElecInS = find(IVC.Intact.hgamma(:,2)<p & IVC.Intact_Scamble.hgamma(:,2)<p);
    %         respElecSnI = find(IVC.Scamble.hgamma(:,2)<p & IVC.Scamble_Intact.hgamma(:,2)<p);
    %         respElecN = setdiff(1:size(IVC.Intact.hgamma,1),[respElecInS;respElecSnI]);
    
    % load electrodes coordinates
    load([subPath 'FsRecon/brain3D/ft_elec.mat'])
    % project electrode to the left hemisphere
    elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1) = elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1).*-1;
    
    for ie = respElecN
        
        hN = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[20 202 18]./255);
    end
    
    for ie = respElecInS'
        
        hInS = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[255 0 0]./255,'fill');
    end
    
    for ie = respElecSnI'
        
        hSnI = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[0 0 255]./255,'fill');
    end
    
    
    
end

title(['BetaIVCdiff(p<' num2str(p) ')'])
legend([hInS,hSnI,hN],['Intact > Scrambled&0'],['Scrambled > Intact&0'],['Other'],'Location','northeast')
saveas(stereoview,[resultPath 'IVCdistribution/' stereoaxes.Title.String '.fig'])
saveas(stereoview,[resultPath 'IVCdistribution/' stereoaxes.Title.String '.jpeg'])

%% plot Power preference electrodes on MNI brain

dirname = uigetdir([resultPath '/avgpow/']);

PlotMNIBrain
allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub


for isub = 1:numel(allsub)%[1,3,4,6,7,8,9,10]%1:numel(allsub)
    subname = allsub{isub};
    subPath = [basePath subname filesep];
    dataPath = [subPath filesep 'Analysis' filesep];
    fprintf(['\n Currently calculating subject: ' subname])
    
    
    p=0.005; % threshold for difference
    %         if exist([dataPath subname 'IVC.mat'],'file')
    %             load([dataPath subname 'IVC']);
    %         end
    
    load([dirname filesep subname])
    
    
    respElecInS = find(pIvS<p & mean(metricM,1)>mean(metricS,1));
    respElecSnI = find(pIvS<p & mean(metricM,1)<mean(metricS,1));
    respElecN = setdiff(1:numel(pIvS),[respElecInS,respElecSnI]);
    
    
    % load electrodes coordinates
    load([subPath 'FsRecon/brain3D/ft_elec.mat'])
    % project electrode to the left hemisphere
    elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1) = elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1).*-1;
    
    for ie = respElecN
        
        hN = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),30,[20 202 18]./255);
    end
    
    for ie = respElecInS'
        
        hInS = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[255 0 0]./255,'fill');
    end
    
    for ie = respElecSnI'
        
        hSnI = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[0 0 255]./255,'fill');
    end
    
    
    
end
title(['ThetaPowDiff(p<' num2str(p) ')'])
legend(stereoaxes,[hInS,hSnI,hN],['Intact > Scrambled'],['Scrambled > Intact'],['Other'],'Location','northeast')
legend('boxoff')
saveas(stereoview,[dirname filesep stereoaxes.Title.String '.fig'])
saveas(stereoview,[dirname filesep stereoaxes.Title.String '.jpeg'])


%% plot Power(low/high freq component) preference electrodes on MNI brain


    allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
        'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub
    
for p=[0.005,0.0025] % threshold for power difference
    
    PlotMNIBrain

    for isub = 1:numel(allsub)%[1,3,4,6,7,8,9,10]%1:numel(allsub)
        subname = allsub{isub};
        subPath = [basePath subname filesep];
        dataPath = [subPath filesep 'Analysis' filesep];
        fprintf(['Currently calculating subject: ' subname '\n'])
        
        
        
        if exist([dataPath subname 'PowDiff.mat'],'file')
            load([dataPath subname 'PowDiff']);
        end

                %%%%---------- original component ----------%%%%
        
%                                 respElecInS = find(PowDiff.theta(:,2)<p & PowDiff.theta(:,1)>0);
%                         respElecSnI = find(PowDiff.theta(:,2)>1-p & PowDiff.theta(:,1)<0);
%                         respElecN = setdiff(1:size(PowDiff.theta,1),[respElecInS;respElecSnI]);
%         
%                                 respElecInS = find(PowDiff.alpha(:,2)<p & PowDiff.alpha(:,1)>0);
%                         respElecSnI = find(PowDiff.alpha(:,2)>1-p & PowDiff.alpha(:,1)<0);
%                         respElecN = setdiff(1:size(PowDiff.alpha,1),[respElecInS;respElecSnI]);
                        
%                                 respElecInS = find(PowDiff.lbeta(:,2)<p & PowDiff.lbeta(:,1)>0);
%                         respElecSnI = find(PowDiff.lbeta(:,2)>1-p & PowDiff.lbeta(:,1)<0);
%                         respElecN = setdiff(1:size(PowDiff.lbeta,1),[respElecInS;respElecSnI]);
%                                 respElecInS = find(PowDiff.hbeta(:,2)<p & PowDiff.hbeta(:,1)>0);
%                         respElecSnI = find(PowDiff.hbeta(:,2)>1-p & PowDiff.hbeta(:,1)<0);
%                         respElecN = setdiff(1:size(PowDiff.hbeta,1),[respElecInS;respElecSnI]);
                        
%                                 respElecInS = find(PowDiff.lgamma(:,2)<p & PowDiff.lgamma(:,1)>0);
%                         respElecSnI = find(PowDiff.lgamma(:,2)>1-p & PowDiff.lgamma(:,1)<0);
%                         respElecN = setdiff(1:size(PowDiff.lgamma,1),[respElecInS;respElecSnI]);
%                                 respElecInS = find(PowDiff.hgamma(:,2)<p & PowDiff.hgamma(:,1)>0);
%                         respElecSnI = find(PowDiff.hgamma(:,2)>1-p & PowDiff.hgamma(:,1)<0);
%                         respElecN = setdiff(1:size(PowDiff.hgamma,1),[respElecInS;respElecSnI]);


        %%%%---------- low frequency component ----------%%%%
        
%                                 respElecInS = find(lpPowDiff.theta(:,2)<p & lpPowDiff.theta(:,1)>0);
%                         respElecSnI = find(lpPowDiff.theta(:,2)>1-p & lpPowDiff.theta(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.theta,1),[respElecInS;respElecSnI]);
%         
%                                 respElecInS = find(lpPowDiff.alpha(:,2)<p & lpPowDiff.alpha(:,1)>0);
%                         respElecSnI = find(lpPowDiff.alpha(:,2)>1-p & lpPowDiff.alpha(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.alpha,1),[respElecInS;respElecSnI]);
%                         
%                                 respElecInS = find(lpPowDiff.lbeta(:,2)<p & lpPowDiff.lbeta(:,1)>0);
%                         respElecSnI = find(lpPowDiff.lbeta(:,2)>1-p & lpPowDiff.lbeta(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.lbeta,1),[respElecInS;respElecSnI]);
%                                 respElecInS = find(lpPowDiff.hbeta(:,2)<p & lpPowDiff.hbeta(:,1)>0);
%                         respElecSnI = find(lpPowDiff.hbeta(:,2)>1-p & lpPowDiff.hbeta(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.hbeta,1),[respElecInS;respElecSnI]);
%                         
%                                 respElecInS = find(lpPowDiff.lgamma(:,2)<p & lpPowDiff.lgamma(:,1)>0);
%                         respElecSnI = find(lpPowDiff.lgamma(:,2)>1-p & lpPowDiff.lgamma(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.lgamma,1),[respElecInS;respElecSnI]);
%                                 respElecInS = find(lpPowDiff.hgamma(:,2)<p & lpPowDiff.hgamma(:,1)>0);
%                         respElecSnI = find(lpPowDiff.hgamma(:,2)>1-p & lpPowDiff.hgamma(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.hgamma,1),[respElecInS;respElecSnI]);
%                
        %%%%---------- High frequency component ----------%%%%
        
%                                 respElecInS = find(lpPowDiff.theta(:,2)<p & lpPowDiff.theta(:,1)>0);
%                         respElecSnI = find(lpPowDiff.theta(:,2)>1-p & lpPowDiff.theta(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.theta,1),[respElecInS;respElecSnI]);
%         
%                                 respElecInS = find(lpPowDiff.alpha(:,2)<p & lpPowDiff.alpha(:,1)>0);
%                         respElecSnI = find(lpPowDiff.alpha(:,2)>1-p & lpPowDiff.alpha(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.alpha,1),[respElecInS;respElecSnI]);
%                         
%                                 respElecInS = find(lpPowDiff.lbeta(:,2)<p & lpPowDiff.lbeta(:,1)>0);
%                         respElecSnI = find(lpPowDiff.lbeta(:,2)>1-p & lpPowDiff.lbeta(:,1)<0);
%                         respElecN = setdiff(1:size(lpPowDiff.lbeta,1),[respElecInS;respElecSnI]);
                                respElecInS = find(hpPowDiff.hbeta(:,2)<p & hpPowDiff.hbeta(:,1)>0);
                        respElecSnI = find(hpPowDiff.hbeta(:,2)>1-p & hpPowDiff.hbeta(:,1)<0);
                        respElecN = setdiff(1:size(hpPowDiff.hbeta,1),[respElecInS;respElecSnI]);
                        
%                                 respElecInS = find(hpPowDiff.lgamma(:,2)<p & hpPowDiff.lgamma(:,1)>0);
%                         respElecSnI = find(hpPowDiff.lgamma(:,2)>1-p & hpPowDiff.lgamma(:,1)<0);
%                         respElecN = setdiff(1:size(hpPowDiff.lgamma,1),[respElecInS;respElecSnI]);
%                                 respElecInS = find(hpPowDiff.hgamma(:,2)<p & hpPowDiff.hgamma(:,1)>0);
%                         respElecSnI = find(hpPowDiff.hgamma(:,2)>1-p & hpPowDiff.hgamma(:,1)<0);
%                         respElecN = setdiff(1:size(hpPowDiff.hgamma,1),[respElecInS;respElecSnI]);
               
        % load electrodes coordinates
        load([subPath 'FsRecon/brain3D/ft_elec.mat'])
        % project electrode to the left hemisphere
        elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1) = elec.MNICoordinate(elec.MNICoordinate(:,1)>0,1).*-1;
        
        for ie = respElecN
            
            hN = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[20 202 18]./255);
        end
        
        for ie = respElecInS'
            
            hInS = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[255 0 0]./255,'fill');
        end
        
        for ie = respElecSnI'
            
            hSnI = scatter3(stereoaxes,elec.MNICoordinate(ie,1),elec.MNICoordinate(ie,2),elec.MNICoordinate(ie,3),60,[0 0 255]./255,'fill');
        end
        
        
        
    end
    
    title(['HighpassLowgammaPowdiff(p<' num2str(2*p) ')'])
    legend([hInS,hSnI,hN],['Intact > Scrambled'],['Scrambled > Intact'],['Other'],'Location','northeast')
    saveas(stereoview,[resultPath 'PowDiff/' stereoaxes.Title.String '.fig'])
    saveas(stereoview,[resultPath 'PowDiff/' stereoaxes.Title.String '.jpeg'])
    close all
    
end
%% plot result electrodes on MNI brain

% PlotMNIBrain


allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub


[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/TFR/*.mat']);

if ~filterindex
    return
end
% load data containing eletrode coordinates
load([pathname filename])

subColor = distinguishable_colors(2*numel(allsub));
subColor = subColor([1 2 3 5 7 9 11 15 16 18],:);
% load in electrodes coordinates
lengendName = {};
nsub = 1;
nelec = 0;
he = [];
for isub = 1:numel(allsub)

    MNICoordinate = importdata([dataPath allsub{isub} '/FsRecon/brain3D/MNI152_coordinates_ras.txt']);
    MNICoordinate = MNICoordinate(:,1:3);
    % choose electrode in ROI
    MNICoordinate = intersect(MNICoordinate,Para.MNICoordinate,'rows');
    MNICoordinate(MNICoordinate(:,1)>0,1) = MNICoordinate(MNICoordinate(:,1)>0,1).*-1;
    if ~isempty(MNICoordinate)
        for ie = 1:size(MNICoordinate,1)
            
            he(nsub) = scatter3(stereoaxes,MNICoordinate(ie,1),MNICoordinate(ie,2),MNICoordinate(ie,3),60,subColor(isub,:),'fill');
        end
        
        lengendName(nsub) = allsub(isub);
        nsub = nsub +1;
        nelec = nelec + size(MNICoordinate,1);
    end
end

 legend(he,lengendName,'Location','northeast')



%% plot all electrodes on MNI brain

PlotMNIBrain


allsub = {'Patient1','Patient2','Patient3','Patient4','Patient6', ...
    'Patient8','Patient9','Patient11','Patient12','Patient13'}; % all sub


subColor = distinguishable_colors(2*numel(allsub));
subColor = subColor([1 2 3 5 7 9 11 15 16 18],:);
% load in electrodes coordinates
lengendName = {};
nsub = 1;
nelec = 0;
he = [];
for isub = 1:numel(allsub)

    MNICoordinate = importdata([dataPath allsub{isub} '/FsRecon/brain3D/MNI152_coordinates_ras.txt']);
    MNICoordinate = MNICoordinate(:,1:3);
    % choose electrode in ROI
    MNICoordinate(MNICoordinate(:,1)>0,1) = MNICoordinate(MNICoordinate(:,1)>0,1).*-1;
    if ~isempty(MNICoordinate)
        for ie = 1:size(MNICoordinate,1)
            
            he(nsub) = scatter3(stereoaxes,MNICoordinate(ie,1),MNICoordinate(ie,2),MNICoordinate(ie,3),60,subColor(isub,:),'fill');
        end
        
        lengendName(nsub) = allsub(isub);
        nsub = nsub +1;
        nelec = nelec + size(MNICoordinate,1);
    end
end

 legend(he,lengendName,'Location','northeast')

%% plot all electrodes on MNI brain （SEEG in Dutch patient）

dataPath = ['~/Desktop/'];

PlotMNIBrain

allsub = {'A','B','C','D','E','F','G'}; % all sub excluding P27


% subColor = distinguishable_colors(numel(allsub));
subColor = [204, 2, 2;...
    85, 147, 3;...
    6, 6, 204;...
    5, 175, 221;...
     112, 42, 3;...
     226, 220, 5;...
     255, 128, 0]./255;
 
% load in electrodes coordinates
lengendName = {};
nsub = 1;
nelec = 0;
he = [];

% Load T1 mask
% aparc_nii = load_nifti([basePath 'Atlas' filesep 'T1excInsula.nii']);
% aparc_vol = round(aparc_nii.vol);
% allroi_index = [];
% temp_indx = find(aparc_vol==1);
% allroi_index = [allroi_index;temp_indx];
% % atlas parcellation coordinates
% [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
% ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
% ras_coordiantes = ras_coordiantes(1:3,:)';
% T1coordiantes = ras_coordiantes;

% Load anatomy atlas
aparc_nii = load_nifti([basePath 'Atlas' filesep 'fsHammers_mith.nii']);
anatomy_vol = round(aparc_nii.vol);
load([basePath 'Atlas' filesep 'Hammers_mith.mat']);
Labels = [{'Unknown'};Labels];


        
fid=fopen([MNI152path '/mri/tempcoordiantes.txt'],'w+');

    allRegion = {};
for isub = 1:numel(allsub)

    load([dataPath allsub{isub} '/brain3D/elec.mat']);
    % choose electrode in ROI
%     tempdev = pdist2(MNICoordinate,T1coordiantes);
%     [it,~] = find(tempdev <=1);
%     MNICoordinate = MNICoordinate(unique(it),:);
    
    VOLCoordinate = aparc_nii.vox2ras\[MNICoordinate ones(size(MNICoordinate,1),1)]';
    VOLCoordinate = round(VOLCoordinate(1:3,:))';
    
    % Project all electrodes to the left hemisphere
%         MNICoordinate(MNICoordinate(:,1)>0,1) = MNICoordinate(MNICoordinate(:,1)>0,1).*-1;
    
    if ~isempty(MNICoordinate)
        for ie = 1:size(MNICoordinate,1)
             tempInd = anatomy_vol(VOLCoordinate(ie,1)-3:VOLCoordinate(ie,1)+3,VOLCoordinate(ie,2)-3:VOLCoordinate(ie,2)+3,...
                 VOLCoordinate(ie,3)-3:VOLCoordinate(ie,3)+3);
           tempInd(tempInd==0) = nan;
           if any(tempInd(:))
             elecRegion = Labels(1+mode(tempInd(:)));
           else
               elecRegion = Labels(1);
           end
           allRegion = [allRegion;elecRegion];
            fprintf(fid,'%s\t%f\t%f\t%f\t%s\n',allsub{isub},MNICoordinate(ie,:),elecRegion{1});
            % plot electrode as 3D sphere
            d = 1;
            [xx,yy,zz] = sphere(20);
            fvc = surf2patch(d*xx+MNICoordinate(ie,1),d*yy+MNICoordinate(ie,2),d*zz+MNICoordinate(ie,3));
            he(nsub) = patch(fvc,'FaceColor',subColor(isub,:),'FaceAlpha',1,'EdgeAlpha', 0);
            
            % plot electrode as 3D scatter dot
%             he(nsub) = scatter3(stereoaxes,MNICoordinate(ie,1),MNICoordinate(ie,2),MNICoordinate(ie,3),60,subColor(isub,:),'fill');
        end
        
        lengendName(nsub) = allsub(isub);
        nsub = nsub +1;
        nelec = nelec + size(MNICoordinate,1);
    end
end

fclose(fid);

gen_roi2warp_code = ['. ~/.zshrc;' ...
    '3dUndump -master ' MNI152path '/mri/T1.nii -orient LPI ' ...
    '-srad 2 -xyz -prefix ' MNI152path '/mri/AllelecMNI152.nii ' ...
    MNI152path '/mri/tempcoordiantes.txt'];
unix(gen_roi2warp_code)
                
 legend(he,lengendName,'Location','northeast')

 allRegionLR = allRegion;
for ir = 1:numel(allRegion)
    if strcmp(allRegion{ir}(end-1:end),' L') | strcmp(allRegion{ir}(end-1:end),' R')
    allRegion{ir}(end-1:end) = [];
    end
end
uniRegion = unique(allRegion);
regionCount = zeros(numel(uniRegion),2);

for ir = 1:numel(allRegionLR)
    if any(strcmp(allRegionLR{ir},uniRegion))
        regionCount(strcmp(allRegionLR{ir},uniRegion),1) = regionCount(strcmp(allRegionLR{ir},uniRegion),1)+1;
        regionCount(strcmp(allRegionLR{ir},uniRegion),2) = nan;
    else
        if any(strcmp(allRegionLR{ir}(1:end-2),uniRegion)) & strcmp(allRegionLR{ir}(end-1:end),' L')
                    regionCount(strcmp(allRegionLR{ir}(1:end-2),uniRegion),1) = regionCount(strcmp(allRegionLR{ir}(1:end-2),uniRegion),1)+1;
        elseif any(strcmp(allRegionLR{ir}(1:end-2),uniRegion)) & strcmp(allRegionLR{ir}(end-1:end),' R')
            regionCount(strcmp(allRegionLR{ir}(1:end-2),uniRegion),2) = regionCount(strcmp(allRegionLR{ir}(1:end-2),uniRegion),2)+1;
        end
    end
end

%% plot specifc regions

 % precentral
aprac_indx = [1];
    bashcode = ['. ~/.zshrc;' ...
        'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
        ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
% bashcode = ['. ~/.zshrc;' ...
%     'fslmaths ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%     ' -thr ' num2str(min(aprac_indx)) ' -uthr ' num2str(max(aprac_indx)) ' ' [MNI152path '/mri/surf_aprac.nii']];
unix(bashcode);

make_outer_surface_wlab ([MNI152path '/mri/surf_aprac.nii'], 15, [MNI152path '/mri/surf_aprac.surf']);

fprintf('Adding surface\n')
[vertices, faces]=read_surf([MNI152path '/mri/surf_aprac.surf']);
faces = faces+1;
aparc_stereo1 = patch(struct(...
    'vertices', vertices, 'faces', faces), ...
    'Parent',stereoaxes, ...
    'FaceColor',[209,73,78]./255, ...
    'FaceAlpha',0, ...
    'EdgeColor', 'none', ...
    'Tag', num2str(aprac_indx));
material dull

% % --------------------------------
 % BA44
% aprac_indx = [1 255];
% bashcode = ['. ~/.zshrc;' ...
%     'mri_binarize --i ' [MNI152path '/mri/fsLeft_Broca_44.nii'] ...
%     ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
% % bashcode = ['. ~/.zshrc;' ...
% %     'fslmaths ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
% %     ' -thr ' num2str(min(aprac_indx)) ' -uthr ' num2str(max(aprac_indx)) ' ' [MNI152path '/mri/surf_aprac.nii']];
% unix(bashcode);
% 
% make_outer_surface_wlab ([MNI152path '/mri/surf_aprac.nii'], 15, [MNI152path '/mri/surf_aprac.surf']);
% 
% fprintf('Adding surface\n')
% [vertices, faces]=read_surf([MNI152path '/mri/surf_aprac.surf']);
% faces = faces+1;
% aparc_stereo2 = patch(struct(...
%     'vertices', vertices, 'faces', faces), ...
%     'Parent',stereoaxes, ...
%     'FaceColor',[255 218 185]./255, ...
%     'FaceAlpha',0.3, ...
%     'EdgeColor', 'none', ...
%     'Tag', num2str(aprac_indx));
% material dull

% % --------------------------------
% supramarginal
aprac_indx = [63];
    bashcode = ['. ~/.zshrc;' ...
        'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
        ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
% bashcode = ['. ~/.zshrc;' ...
%     'fslmaths ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%     ' -thr ' num2str(min(aprac_indx)) ' -uthr ' num2str(max(aprac_indx)) ' ' [MNI152path '/mri/surf_aprac.nii']];
unix(bashcode);

make_outer_surface_wlab ([MNI152path '/mri/surf_aprac.nii'], 15, [MNI152path '/mri/surf_aprac.surf']);

fprintf('Adding surface\n')
[vertices, faces]=read_surf([MNI152path '/mri/surf_aprac.surf']);
faces = faces+1;
aparc_stereo3 = patch(struct(...
    'vertices', vertices, 'faces', faces), ...
    'Parent',stereoaxes, ...
    'FaceColor',[238 130 238]./255, ...
    'FaceAlpha',0, ...
    'EdgeColor', 'none', ...
    'Tag', num2str(aprac_indx));
material dull

% --------------------------
% middle occipital
aprac_indx = [51];
    bashcode = ['. ~/.zshrc;' ...
        'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
        ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
% bashcode = ['. ~/.zshrc;' ...
%     'fslmaths ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%     ' -thr ' num2str(min(aprac_indx)) ' -uthr ' num2str(max(aprac_indx)) ' ' [MNI152path '/mri/surf_aprac.nii']];
unix(bashcode);

make_outer_surface_wlab ([MNI152path '/mri/surf_aprac.nii'], 15, [MNI152path '/mri/surf_aprac.surf']);

fprintf('Adding surface\n')
[vertices, faces]=read_surf([MNI152path '/mri/surf_aprac.surf']);
faces = faces+1;
aparc_stereo4 = patch(struct(...
    'vertices', vertices, 'faces', faces), ...
    'Parent',stereoaxes, ...
    'FaceColor',[18 53 85]./255, ...
    'FaceAlpha',0, ...
    'EdgeColor', 'none', ...
    'Tag', num2str(aprac_indx));
material dull



% legend([aparc_stereo1,aparc_stereo3,aparc_stereo4],['Precentral (8sub,59 elec)'],...
% ['SupraMarginal (9sub,56 elec)'],['MiddleOccipital (6 sub,40 elec)'],'Location','northeast')

% abbreviation version
legend([aparc_stereo1,aparc_stereo3,aparc_stereo4],['PreCG (8sub,59 elec)'],...
['SMG (9sub,56 elec)'],['MOG (6 sub,40 elec)'],'Location','northeast')

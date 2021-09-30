%% data path definition
basePath = '~/Desktop/ActionPrediction/';
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
    elec.elecposMNI(elec.elecposMNI(:,1)>0,1) = elec.elecposMNI(elec.elecposMNI(:,1)>0,1).*-1;
    
    for ie = respElecI'
        
        hI(ie) = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[255 0 0]./255,'fill');
    end
    
    for ie = respElecS'
        
        hS(ie) = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[0 0 255]./255,'fill');
    end
    
    for ie = respElecInS'
        
        hInS(ie) = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[0 255 0]./255,'fill');
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
    elec.elecposMNI(elec.elecposMNI(:,1)>0,1) = elec.elecposMNI(elec.elecposMNI(:,1)>0,1).*-1;
    
    for ie = respElecN
        
        hN = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[20 202 18]./255);
    end
    
    for ie = respElecInS'
        
        hInS = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[255 0 0]./255,'fill');
    end
    
    for ie = respElecSnI'
        
        hSnI = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[0 0 255]./255,'fill');
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
    elec.elecposMNI(elec.elecposMNI(:,1)>0,1) = elec.elecposMNI(elec.elecposMNI(:,1)>0,1).*-1;
    
    for ie = respElecN
        
        hN = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),30,[20 202 18]./255);
    end
    
    for ie = respElecInS'
        
        hInS = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[255 0 0]./255,'fill');
    end
    
    for ie = respElecSnI'
        
        hSnI = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[0 0 255]./255,'fill');
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
        elec.elecposMNI(elec.elecposMNI(:,1)>0,1) = elec.elecposMNI(elec.elecposMNI(:,1)>0,1).*-1;
        
        for ie = respElecN
            
            hN = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[20 202 18]./255);
        end
        
        for ie = respElecInS'
            
            hInS = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[255 0 0]./255,'fill');
        end
        
        for ie = respElecSnI'
            
            hSnI = scatter3(stereoaxes,elec.elecposMNI(ie,1),elec.elecposMNI(ie,2),elec.elecposMNI(ie,3),60,[0 0 255]./255,'fill');
        end
        
        
        
    end
    
    title(['HighpassLowgammaPowdiff(p<' num2str(2*p) ')'])
    legend([hInS,hSnI,hN],['Intact > Scrambled'],['Scrambled > Intact'],['Other'],'Location','northeast')
    saveas(stereoview,[resultPath 'PowDiff/' stereoaxes.Title.String '.fig'])
    saveas(stereoview,[resultPath 'PowDiff/' stereoaxes.Title.String '.jpeg'])
    close all
    
end
%% plot result electrodes on MNI brain

[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/TFR/*.mat']);

if ~filterindex
    return
end

% load in TFR LMEM
load([pathname filename])
Para.elecposMNI(Para.elecposMNI(:,1)>0,1) = Para.elecposMNI(Para.elecposMNI(:,1)>0,1).*-1;
for ie = 1:size(Para.elecposMNI,1)
    
    he(ie) = scatter3(stereoaxes,Para.elecposMNI(ie,1),Para.elecposMNI(ie,2),Para.elecposMNI(ie,3),60,[255 185 15]./255,'fill');
end


% 240,0,0 precentral
% 193 255 193 Pft
% 100 149 237 Visual
% 46 139 87 BA44
% 255 185 15 IPL
%% plot specifc regions

%  % precentral
% aprac_indx = [1];
%     bashcode = ['. ~/.zshrc;' ...
%         'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%         ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
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
% aparc_stereo = patch(struct(...
%     'vertices', vertices, 'faces', faces), ...
%     'Parent',stereoaxes, ...
%     'FaceColor',[255 218 185]./255, ...
%     'FaceAlpha',0.5, ...
%     'EdgeColor', 'none', ...
%     'Tag', num2str(aprac_indx));
% material dull
%
% --------------------------
% aprac_indx = [51];
%     bashcode = ['. ~/.zshrc;' ...
%         'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%         ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
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
% aparc_stereo = patch(struct(...
%     'vertices', vertices, 'faces', faces), ...
%     'Parent',stereoaxes, ...
%     'FaceColor',[144 238 144]./255, ...
%     'FaceAlpha',0.5, ...
%     'EdgeColor', 'none', ...
%     'Tag', num2str(aprac_indx));
% material dull
%
% aprac_indx = [53];
%     bashcode = ['. ~/.zshrc;' ...
%         'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%         ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
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
% aparc_stereo = patch(struct(...
%     'vertices', vertices, 'faces', faces), ...
%     'Parent',stereoaxes, ...
%     'FaceColor',[30 144 255]./255, ...
%     'FaceAlpha',0.5, ...
%     'EdgeColor', 'none', ...
%     'Tag', num2str(aprac_indx));
% material dull
%
%
% aprac_indx = [57];
%     bashcode = ['. ~/.zshrc;' ...
%         'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%         ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
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
% aparc_stereo = patch(struct(...
%     'vertices', vertices, 'faces', faces), ...
%     'Parent',stereoaxes, ...
%     'FaceColor',[255 246 143]./255, ...
%     'FaceAlpha',0.5, ...
%     'EdgeColor', 'none', ...
%     'Tag', num2str(aprac_indx));
% material dull
%
%
% aprac_indx = [61];
%     bashcode = ['. ~/.zshrc;' ...
%         'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%         ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
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
% aparc_stereo = patch(struct(...
%     'vertices', vertices, 'faces', faces), ...
%     'Parent',stereoaxes, ...
%     'FaceColor',[255 106 106]./255, ...
%     'FaceAlpha',0.5, ...
%     'EdgeColor', 'none', ...
%     'Tag', num2str(aprac_indx));
% material dull
%
%
% % another region
% aprac_indx = [63];
%     bashcode = ['. ~/.zshrc;' ...
%         'mri_binarize --i ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%         ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
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
% aparc_stereo = patch(struct(...
%     'vertices', vertices, 'faces', faces), ...
%     'Parent',stereoaxes, ...
%     'FaceColor',[238 130 238]./255, ...
%     'FaceAlpha',0.5, ...
%     'EdgeColor', 'none', ...
%     'Tag', num2str(aprac_indx));
% material dull
% % --------------------------------


%  % pft
aprac_indx = [1 255];
bashcode = ['. ~/.zshrc;' ...
    'mri_binarize --i ' [MNI152path '/mri/fsLeft_IPL_PFt.nii'] ...
    ' --min ' num2str(min(aprac_indx)-0.01) ' --max ' num2str(max(aprac_indx)+0.01) ' --o ' [MNI152path '/mri/surf_aprac.nii']];
% bashcode = ['. ~/.zshrc;' ...
%     'fslmaths ' [MNI152path '/mri/fsAnatomyMacro.nii'] ...
%     ' -thr ' num2str(min(aprac_indx)) ' -uthr ' num2str(max(aprac_indx)) ' ' [MNI152path '/mri/surf_aprac.nii']];
unix(bashcode);

make_outer_surface_wlab ([MNI152path '/mri/surf_aprac.nii'], 15, [MNI152path '/mri/surf_aprac.surf']);

fprintf('Adding surface\n')
[vertices, faces]=read_surf([MNI152path '/mri/surf_aprac.surf']);
faces = faces+1;
aparc_stereo = patch(struct(...
    'vertices', vertices, 'faces', faces), ...
    'Parent',stereoaxes, ...
    'FaceColor',[255 218 185]./255, ...
    'FaceAlpha',0.5, ...
    'EdgeColor', 'none', ...
    'Tag', num2str(aprac_indx));
material dull

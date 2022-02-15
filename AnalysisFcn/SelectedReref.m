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


%% -------- ROI Partition --------

% Region of Interest include:
ROIIndex = {[1,2],[49,50],[51,52],[53,54],[59,60],[61,62],[63,64],[1],[1]};
ROIAtlas = {'fsAnatomyMacro.nii','fsAnatomyMacro.nii','fsAnatomyMacro.nii','fsAnatomyMacro.nii', ...
    'fsAnatomyMacro.nii','fsAnatomyMacro.nii','fsAnatomyMacro.nii','BA44.nii','PFt.nii'};
ROIText = {'Precentral','SuperiorOccipitalGyrus','MiddleOccipitalGyrus',...
    'InferiorOccipitalGyrus','SuperiorParietalLobe','InferiorParietalLobe','SupraMarginal','BA44','PFt'};

roiDist = 1; % maximum distance between electrodes and ROI voxels

seedIndex = [1];
searchIndex = [7];

 iPreCG = 1;
 iSMG = 7;
 iMOG = 3;
        
        % load altlas infomation for seed roi
        aparc_nii = load_nifti([basePath 'Atlas' filesep ROIAtlas{iPreCG}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        if ndims(aparc_vol) == 3
            for iroi = ROIIndex{iPreCG}
                temp_indx = find(aparc_vol==iroi);
                allroi_index = [allroi_index;temp_indx];
            end
        elseif ndims(aparc_vol) == 4
            %%%------ For probabilistic map ------%%%
            aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{iPreCG}),4));
            allroi_index = find(aparc_vol>=10);  % minimum probability
        end
        
        % atlas parcellation coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        PreCG_coordiantes = ras_coordiantes;
        
        % load altlas infomation for search roi
        aparc_nii = load_nifti([basePath 'Atlas' filesep ROIAtlas{iSMG}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        if ndims(aparc_vol) == 3
            for iroi = ROIIndex{iSMG}
                temp_indx = find(aparc_vol==iroi);
                allroi_index = [allroi_index;temp_indx];
            end
        elseif ndims(aparc_vol) == 4
            %%%------ For probabilistic map ------%%%
            aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{iSMG}),4));
            allroi_index = find(aparc_vol>=10);  % minimum probability
        end
        
        % atlas parcellation coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        SMG_coordiantes = ras_coordiantes;

        
         % load altlas infomation for search roi
        aparc_nii = load_nifti([basePath 'Atlas' filesep ROIAtlas{iMOG}]);
        aparc_vol = round(aparc_nii.vol);
        allroi_index = [];
        %%%------ For maxprobabilistic map ------%%%
        if ndims(aparc_vol) == 3
            for iroi = ROIIndex{iMOG}
                temp_indx = find(aparc_vol==iroi);
                allroi_index = [allroi_index;temp_indx];
            end
        elseif ndims(aparc_vol) == 4
            %%%------ For probabilistic map ------%%%
            aparc_vol = squeeze(sum(aparc_vol(:,:,:,ROIIndex{iMOG}),4));
            allroi_index = find(aparc_vol>=10);  % minimum probability
        end
        
        % atlas parcellation coordinates
        [tx,ty,tj] = ind2sub(size(aparc_vol),allroi_index);
        ras_coordiantes = aparc_nii.vox2ras*[[tx,ty,tj] ones(size(tx,1),1)]';
        ras_coordiantes = ras_coordiantes(1:3,:)';
        MOG_coordiantes = ras_coordiantes;      
        
        
        %% -------- Subject Level --------
        

            subname = patientID;
            subPath = [basePath filesep 'Data' filesep subname filesep];
            dataPath = [subPath filesep 'Analysis' filesep];
            fprintf(['\n Currently calculating subject: ' subname])
            
            % choose seed electrodes according to MNI coordinates
            elecposMNI = ftData.elec.elecposMNI;
            tempdev = pdist2(elecposMNI,PreCG_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            PreCGElec = unique(it);
            % choose search electrodes according to MNI coordinates
            tempdev = pdist2(elecposMNI,SMG_coordiantes);
            [it,~] = find(tempdev <=roiDist);
           SMGElec = unique(it);
            % choose search electrodes according to MNI coordinates
            tempdev = pdist2(elecposMNI,MOG_coordiantes);
            [it,~] = find(tempdev <=roiDist);
            MOGElec = unique(it);

            
            for i = PreCGElec'
                for j = SMGElec'
                    overlap = find(rerefMap(i,:)~=0&rerefMap(j,:)~=0);
                    if ~isempty(overlap)
                        for io = overlap'
                                                        disp(io)
                            if distances(i,io)<=distances(j,io)
                                rerefMap(j,io) = 0;
                            else
                                rerefMap(i,io) = 0;
                            end
                        end
                    end
                end
            end
            
            for i = SMGElec'
                for j = MOGElec'
                    overlap = find(rerefMap(i,:)~=0&rerefMap(j,:)~=0);
                    if ~isempty(overlap)
                        for io = overlap'
                            disp(io)
                            if distances(i,io)<=distances(j,io)
                                rerefMap(j,io) = 0;
                            else
                                rerefMap(i,io) = 0;
                            end
                        end
                    end
                end
            end
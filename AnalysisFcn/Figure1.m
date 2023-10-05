%% User defined data path
dataPath = 'C:\Users\qin2\Documents\ActionPredictionECoG\Data\ECoGData\';

%% plot all electrodes on MNI brain
addpath('Atlas')
addpath('tools')

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
    subname = allsub{isub};
    load([dataPath subname 'LARER_rerefDataInROI.mat'])
    MNICoordinate = rerefData.elec.elecposMNI;
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
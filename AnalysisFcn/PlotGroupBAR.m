%% plot ERP
clear;
dirname = uigetdir(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/BPall']);

% load in TFR LMEM
filenames = dir([dirname filesep '*mat']);

barH = [];
barE = [];
barX = [];
xtick = [];
xticklabel = [];
pMatrix = nan(2*size(filenames,1));
for n = 1:size(filenames,1)
    
    load([filenames(n).folder filesep filenames(n).name],'y2plot','se2plot','pMap')
    barH = [barH,y2plot'];
    barE = [barE,se2plot'];
    barX = [barX,[3*(n-1)+1,3*(n-1)+2]];
    xtick = [xtick,mean([3*(n-1)+1,3*(n-1)+2])];
    xticklabel = [xticklabel,{filenames(n).name(1:end-4)}];
    pMatrix(2*n-1,2*n) = pMap;
end

hh = goodbar(barH,barE,pMatrix,barX);
set(gca,'xtick',xtick)
set(gca,'xticklabel',xticklabel)
for ip = 2:2:numel(hh)
    set(hh(ip),'FaceColor',[0.5 0.5 0.5]);
end
legend([hh(1),hh(2)],'Intact','Scrambled')

ylabel('Normalised Power')
% ylim([0.7 1.2])

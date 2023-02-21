%% Quantify the Motion Energy in the Action Prediction Task

clear;clc;
%% Set Directories

movDir = 'C:\Users\qin2\Documents\ActionPredictionECoG\Data\Patient1\Stimuli'; % Needs to be editted. The directory of the folder that contains the videos to be processed.


%% Get the Video List

movListA = readmatrix([movDir filesep 'StimulusList_A.txt'],'OutputType','string'); % Needs to be editted. The format of the videos to be processed.
movListB = readmatrix([movDir filesep 'StimulusList_B.txt'],'OutputType','string'); % Needs to be editted. The format of the videos to be processed.
movListC = readmatrix([movDir filesep 'StimulusList_C.txt'],'OutputType','string'); % Needs to be editted. The format of the videos to be processed.
movA = VideoReader([[movDir filesep 'RunA_squareAtBeginningOnly.avi']]);
movB = VideoReader([[movDir filesep 'RunB_squareAtBeginningOnly.avi']]);
movC = VideoReader([[movDir filesep 'RunC_squareAtBeginningOnly.avi']]);

%% Calculate Motion Energy
timeWin = [-0.5 1];

%% Movie A
strlen = 0;
for imov = 1:numel(movListA)
    tmpStimInfo = readtable([movDir filesep char(movListA(imov)) '.csv']);
    camChgTime = tmpStimInfo.Start_s_(2:end);
    nRows = movA.Height;
    nCols = movA.Width;
    
    % movie wise motion energy
    movStartFrame = round(movA.FrameRate*tmpStimInfo.Start_s_(2))+1;
    movEndFrame = movStartFrame+round(movA.FrameRate*tmpStimInfo.TotalLength(1))-2;
    frames = im2single(read(movA,[movStartFrame movEndFrame]));
    for iframe = 1:size(frames,4)-1
        s = ['Processing Videos ' num2str(imov) '/' num2str(numel(movListA))];
        strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
        strlen = strlentmp - strlen;
        
        energyMov.(movListA(imov,1))(iframe) = mean(vecnorm(reshape(frames(:,:,:,iframe), ...
            nRows*nCols,3)-reshape(frames(:,:,:,iframe+1),nRows*nCols,3),2,2),'omitnan');
    end
  
    % act wise motion energy
    for icam = 1:numel(camChgTime)
        actStartFrame = floor((camChgTime(icam)+timeWin(1))*movA.FrameRate);
        actEndFrame = actStartFrame + round((sum(abs(timeWin)))*movA.FrameRate);
        frames = im2single(read(movA,[actStartFrame actEndFrame]));
        
        for iframe = 1:size(frames,4)-1
            s = ['Processing Video Clips ' num2str(imov) '/' num2str(numel(movListA)) ...
                ',Scene ' num2str(icam) '/' num2str(numel(camChgTime)) ...
                ',Frame ' num2str(iframe) '/' num2str(size(frames,4)-1)];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;

            energy.(movListA(imov,1))(icam,iframe) = mean(vecnorm(reshape(frames(:,:,:,iframe), ...
                nRows*nCols,3)-reshape(frames(:,:,:,iframe+1),nRows*nCols,3),2,2),'omitnan');
        end
    end
    
end


%% Movie B
strlen = 0;
for imov = 1:numel(movListB)
    tmpStimInfo = readtable([movDir filesep char(movListB(imov)) '.csv']);
    camChgTime = tmpStimInfo.Start_s_(2:end);
    nRows = movB.Height;
    nCols = movB.Width;
    
    % movie wise motion energy
    movStartFrame = round(movB.FrameRate*tmpStimInfo.Start_s_(2))+1;
    movEndFrame = movStartFrame+round(movB.FrameRate*tmpStimInfo.TotalLength(1))-2;
    frames = im2single(read(movB,[movStartFrame movEndFrame]));
    for iframe = 1:size(frames,4)-1
        s = ['Processing Videos ' num2str(imov) '/' num2str(numel(movListB))];
        strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
        strlen = strlentmp - strlen;
        
        energyMov.(movListB(imov,1))(iframe) = mean(vecnorm(reshape(frames(:,:,:,iframe), ...
            nRows*nCols,3)-reshape(frames(:,:,:,iframe+1),nRows*nCols,3),2,2),'omitnan');
    end
  
    % act wise motion energy
    for icam = 1:numel(camChgTime)
        actStartFrame = floor((camChgTime(icam)+timeWin(1))*movB.FrameRate);
        actEndFrame = actStartFrame + round((sum(abs(timeWin)))*movB.FrameRate);
        frames = im2single(read(movB,[actStartFrame actEndFrame]));
        
        for iframe = 1:size(frames,4)-1
            s = ['Processing Video Clips ' num2str(imov) '/' num2str(numel(movListB)) ...
                ',Scene ' num2str(icam) '/' num2str(numel(camChgTime)) ...
                ',Frame ' num2str(iframe) '/' num2str(size(frames,4)-1)];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;

            energy.(movListB(imov,1))(icam,iframe) = mean(vecnorm(reshape(frames(:,:,:,iframe), ...
                nRows*nCols,3)-reshape(frames(:,:,:,iframe+1),nRows*nCols,3),2,2),'omitnan');
        end
    end
    
end


%% Movie C
strlen = 0;
for imov = 1:numel(movListC)
    tmpStimInfo = readtable([movDir filesep char(movListC(imov)) '.csv']);
    camChgTime = tmpStimInfo.Start_s_(2:end);
    nRows = movC.Height;
    nCols = movC.Width;
    
    % movie wise motion energy
    movStartFrame = round(movC.FrameRate*tmpStimInfo.Start_s_(2))+1;
    movEndFrame = movStartFrame+round(movC.FrameRate*tmpStimInfo.TotalLength(1))-2;
    frames = im2single(read(movC,[movStartFrame movEndFrame]));
    for iframe = 1:size(frames,4)-1
        s = ['Processing Videos ' num2str(imov) '/' num2str(numel(movListC))];
        strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
        strlen = strlentmp - strlen;
        
        energyMov.(movListC(imov,1))(iframe) = mean(vecnorm(reshape(frames(:,:,:,iframe), ...
            nRows*nCols,3)-reshape(frames(:,:,:,iframe+1),nRows*nCols,3),2,2),'omitnan');
    end
  
    % act wise motion energy
    for icam = 1:numel(camChgTime)
        actStartFrame = floor((camChgTime(icam)+timeWin(1))*movC.FrameRate);
        actEndFrame = actStartFrame + round((sum(abs(timeWin)))*movC.FrameRate);
        frames = im2single(read(movC,[actStartFrame actEndFrame]));
        
        for iframe = 1:size(frames,4)-1
            s = ['Processing Video Clips ' num2str(imov) '/' num2str(numel(movListC)) ...
                ',Scene ' num2str(icam) '/' num2str(numel(camChgTime)) ...
                ',Frame ' num2str(iframe) '/' num2str(size(frames,4)-1)];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;

            energy.(movListC(imov,1))(icam,iframe) = mean(vecnorm(reshape(frames(:,:,:,iframe), ...
                nRows*nCols,3)-reshape(frames(:,:,:,iframe+1),nRows*nCols,3),2,2),'omitnan');
        end
    end
    
end




%% Save
FrameRate = movA.FrameRate;
save([movDir filesep 'MotionEnergy.mat'],'energy','energyMov','timeWin','FrameRate')


%% Analyse motion energy

timePT = min(timeWin)+0.02:1/FrameRate:max(timeWin);
energyM = [];
energyS = [];
energyMovM = [];
energyMovS = [];
movName = fieldnames(energy);
for imov = 1:numel(movName)
    if strncmp(movName{imov},'S',1)
        energyS = [energyS;energy.(movName{imov})];
        energyMovS = [energyMovS;mean(energyMov.(movName{imov}))];
        Mind = strcmp(['M' movName{imov}(2:end)],movName);
        energyM = [energyM;energy.(movName{Mind})];
        energyMovM = [energyMovM;mean(energyMov.(movName{Mind}))];
    end
end

% test difference between conditions
[~,pMapmov,~,tStat] = ttest(energyMovM,energyMovS);
[bf10,p] = bf.ttest2(energyMovM,energyMovS);
[~,pMap,~,~] = ttest2(energyM,energyS);
[p_fdr, p_masked] = fdr(pMap, 0.05);
[FDR] = mafdr(pMap,'BHFDR',true);
% [FDR,Q] = myfdr(pMap,0.05);
highlight = p_masked;
% the significant points could be marked with a red stars
highlight =double(highlight);
highlight(highlight==0) = nan;


hf = figure;
hold on

mean2plot = mean(energyM,1);
ste2plot = std(energyM,0,1)/sqrt(size(energyM,1));
mm = 0;
ss = 0;
tt = min(timePT);
for it = 1:numel(timePT)-1
    mm = [mm,mean2plot(it)];
    ss = [ss,0];
    tt = [tt,timePT(it)+0.001];
    
    mm = [mm,mean2plot(it)];
    ss = [ss,ste2plot(it)];
    tt = [tt,timePT(it)+0.002];
    
    mm = [mm,mean2plot(it)];
    ss = [ss,ste2plot(it)];
    tt = [tt,timePT(it+1)];
    
        mm = [mm,mean2plot(it)];
    ss = [ss,0];
    tt = [tt,timePT(it+1)+0.001];
end

hM = shadedErrorBar(tt, mm,ss,{'color',[255 106 106]/255,'linewidth',1,'linestyle','-'},1);

mean2plot = mean(energyS,1);
ste2plot = std(energyS,0,1)/sqrt(size(energyS,1));
mm = 0;
ss = 0;
tt = min(timePT);
for it = 1:numel(timePT)-1
    mm = [mm,mean2plot(it)];
    ss = [ss,0];
    tt = [tt,timePT(it)+0.001];
    
    mm = [mm,mean2plot(it)];
    ss = [ss,ste2plot(it)];
    tt = [tt,timePT(it)+0.002];
    
    mm = [mm,mean2plot(it)];
    ss = [ss,ste2plot(it)];
    tt = [tt,timePT(it+1)];
    
        mm = [mm,mean2plot(it)];
    ss = [ss,0];
    tt = [tt,timePT(it+1)+0.001];
end

hS = shadedErrorBar(tt, mm, ss,{'color',[30 144 255]/255,'linewidth',1,'linestyle','-'},1);
% hM = bar(timePT+0.02,mean(energyM,1),1,'FaceColor','none','EdgeColor',[255 106 106]/255,'FaceAlpha',1,'EdgeAlpha',0.75,'linewidth',1.5);
% hS = bar(timePT+0.02,mean(energyS,1),1,'FaceColor','none','EdgeColor',[30 144 255]/255,'FaceAlpha',1,'EdgeAlpha',0.75,'linewidth',1.5);
hsig = plot(timePT,1.05*max(mean(energyM,1))*highlight,'k*');
% rectangle('Position',[-0.1 0 0.2 0.4])

xlabel('Time relative to camera change (sec)')
ylabel('Motion Energy')
xlim([-0.48 1])
% plot([0 0],get(gca,'ylim'),'k--')

% legend([hM.mainLine,hS.mainLine,hsig], ...
%     ['Intact'],['Scrambled'],['p<0.05(corrected)'], ...
%     'Location','NorthWest','box','off');
hMsamp = plot([2,3],[1,1],'Color',[255 106 106]/255,'linewidth',1.5);
hSsamp = plot([2,3],[1,1],'Color',[30 144 255]/255,'linewidth',1.5);
legend([hMsamp,hSsamp], ...
    ['Intact'],['Scrambled'], ...
    'Location','NorthWest','box','off');

haZoom = axes(hf,'Position',[0.6 0.6 0.25 0.25],...
    'Xlim',[0.5 2.5],'XTick', [1,2],'XTickLabel',{'I','S'},...
    'Ylim',[0.01 0.045],'FontSize',12,'NextPlot','add');
plot(haZoom,1+zeros(size(energyMovM,1),1),energyMovM,'v','color',[255 106 106]/255)
plot(haZoom,2+zeros(size(energyMovS,1),1),energyMovS,'v','color',[30 144 255]/255)
plot(haZoom,repmat([1.1;1.9],1,numel(energyMovS)),[energyMovM,energyMovS]','k')
text(haZoom,-0.2,1.15,['t_(_1_1_)=' num2str(round(tStat.tstat,2)) ...
    ', BF_1_0=' num2str(round(bf10,2)) ', p=' num2str(round(pMapmov,2))],'unit','Normalized')

saveas(hf,[movDir filesep 'MotionEnergy.fig'])




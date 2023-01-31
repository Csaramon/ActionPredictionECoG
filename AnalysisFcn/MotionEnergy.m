%% Quantify the Motion Energy in the Action Prediction Task


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
movEndFrame = 1;
strlen = 0;
for imov = 1:numel(movListA)
    tmpStimInfo = readtable([movDir filesep char(movListA(imov)) '.csv']);
    camChgTime = tmpStimInfo.Stop_s_(2:end);
    
    movStartFrame = movEndFrame;
    movEndFrame = movStartFrame + round((tmpStimInfo.TotalLength(1)+3) * movA.FrameRate); % Add inter movie interval of 3 seconds
    frames = im2double(read(movA,10*movA.FrameRate+[movStartFrame movEndFrame]));
    nRows = size(frames,1);
    nCols = size(frames,2);
    
    for icam = 1:numel(camChgTime)

        startFrame = round([camChgTime(icam)+timeWin(1)]*movA.FrameRate);
        numFrame = round((sum(abs(timeWin)))*movA.FrameRate);

        for iframe = 1:numFrame
            s = ['Processing Video ' num2str(imov) '/' num2str(numel(movListA)) ...
                ',Scene ' num2str(icam) '/' num2str(numel(camChgTime)) ...
                ',Frame ' num2str(iframe) '/' num2str(numFrame)];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;
%             for row = 1:nRows
%                 for col = 1:nCols
%                     diffMat(row,col) = pdist([reshape(frames(row,col,:,startFrame+iframe-1),1,3);reshape(frames(row,col,:,startFrame+iframe),1,3)],...
%                         'Euclidean');
%                 end
%             end            
                energy.(movListA(imov,1))(icam,iframe) = mean(vecnorm(reshape(frames(:,:,:,startFrame+iframe-1), ...
                    nRows*nCols,3)-reshape(frames(:,:,:,startFrame+iframe),nRows*nCols,3),2,2),'omitnan');
        end
    end

end


%% Movie B
movEndFrame = 1;
strlen = 0;
for imov = 1:numel(movListB)
    tmpStimInfo = readtable([movDir filesep char(movListB(imov)) '.csv']);
    camChgTime = tmpStimInfo.Stop_s_(2:end);
    
    movStartFrame = movEndFrame;
    movEndFrame = movStartFrame + round((tmpStimInfo.TotalLength(1)+3) * movB.FrameRate); % Add inter movie interval of 3 seconds
    frames = im2double(read(movB,10*movB.FrameRate+[movStartFrame movEndFrame]));
    nRows = size(frames,1);
    nCols = size(frames,2);
    
    for icam = 1:numel(camChgTime)

        startFrame = round([camChgTime(icam)+timeWin(1)]*movB.FrameRate);
        numFrame = round((sum(abs(timeWin)))*movB.FrameRate);

        for iframe = 1:numFrame
            s = ['Processing Video ' num2str(imov) '/' num2str(numel(movListB)) ...
                ',Scene ' num2str(icam) '/' num2str(numel(camChgTime)) ...
                ',Frame ' num2str(iframe) '/' num2str(numFrame)];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;
%             for row = 1:nRows
%                 for col = 1:nCols
%                     diffMat(row,col) = pdist([reshape(frames(row,col,:,startFrame+iframe-1),1,3);reshape(frames(row,col,:,startFrame+iframe),1,3)],...
%                         'Euclidean');
%                 end
%             end            
                energy.(movListB(imov,1))(icam,iframe) = mean(vecnorm(reshape(frames(:,:,:,startFrame+iframe-1), ...
                    nRows*nCols,3)-reshape(frames(:,:,:,startFrame+iframe),nRows*nCols,3),2,2),'omitnan');
        end
    end

end


%% Movie C
movEndFrame = 1;
strlen = 0;
for imov = 1:numel(movListC)
    tmpStimInfo = readtable([movDir filesep char(movListC(imov)) '.csv']);
    camChgTime = tmpStimInfo.Stop_s_(2:end);
    
    movStartFrame = movEndFrame;
    movEndFrame = movStartFrame + round((tmpStimInfo.TotalLength(1)+3) * movC.FrameRate); % Add inter movie interval of 3 seconds
    frames = im2double(read(movC,10*movC.FrameRate+[movStartFrame movEndFrame]));
    nRows = size(frames,1);
    nCols = size(frames,2);
    
    for icam = 1:numel(camChgTime)

        startFrame = round([camChgTime(icam)+timeWin(1)]*movC.FrameRate);
        numFrame = round((sum(abs(timeWin)))*movC.FrameRate);

        for iframe = 1:numFrame
            s = ['Processing Video ' num2str(imov) '/' num2str(numel(movListC)) ...
                ',Scene ' num2str(icam) '/' num2str(numel(camChgTime)) ...
                ',Frame ' num2str(iframe) '/' num2str(numFrame)];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;
%             for row = 1:nRows
%                 for col = 1:nCols
%                     diffMat(row,col) = pdist([reshape(frames(row,col,:,startFrame+iframe-1),1,3);reshape(frames(row,col,:,startFrame+iframe),1,3)],...
%                         'Euclidean');
%                 end
%             end            
                energy.(movListC(imov,1))(icam,iframe) = mean(vecnorm(reshape(frames(:,:,:,startFrame+iframe-1), ...
                    nRows*nCols,3)-reshape(frames(:,:,:,startFrame+iframe),nRows*nCols,3),2,2),'omitnan');
        end
    end

end


%% Save
FrameRate = movA.FrameRate;
save([movDir filesep 'MotionEnergy.mat'],'energy','timeWin','FrameRate')


%% Analyse motion energy

timePT = min(timeWin)+0.02:1/FrameRate:max(timeWin);
energyM = [];
energyS = [];
movName = fieldnames(energy);
for imov = 1:numel(movName)
    if strncmp(movName{imov},'S',1)
        energyS = [energyS;energy.(movName{imov})];
    elseif strncmp(movName{imov},'M',1)
        energyM = [energyM;energy.(movName{imov})];
    end
end

% test difference between conditions
[~,pMap] = ttest2(energyM,energyS);
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;
% the significant points could be marked with a red stars
highlight =double(highlight);
highlight(highlight==0) = nan;


figure;
hold on
hM = shadedErrorBar(timePT, mean(energyM,1),std(energyM,0,1)/sqrt(size(energyM,1)),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(timePT, mean(energyS,1), std(energyS,0,1)/sqrt(size(energyS,1)),{'color',[30 144 255]/255},1);
hsig = plot(timePT,1.05*max(mean(energyM,1))*highlight,'k*');

xlabel('Time relative to camera change (sec)')
ylabel('Motion Energy')
xlim([-0.5 1])
plot([0 0],get(gca,'ylim'),'k--')

%% plot ERP
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/ERP/*.mat']);

if ~filterindex 
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end 
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

subplot(2,1,1);
hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
legend([hM.mainLine,hS.mainLine],['Intact'],['Scrambled']);
title(filename(1:end-4))
ylabel('Potential (μV)')
ylim([-5 5])


subplot(2,1,2);
hold on
h = plot(Para.time,tMap);
hsig = plot(Para.time,(min(tMap)-0.1*range(tMap))*highlight,'r*');
text(gca,0.55,3.5, ['Nsub:' num2str(numel(unique(lmeTBL.Sub))) ...
    ' Nelec:' num2str(numel(unique(lmeTBL.Elec)))]);
xlabel('Time relative to camera change (sec)')
ylabel('t Value')
ylim([-4 4])


%% plot Bandpower
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/BP/*.mat']);

if ~filterindex 
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
Ls = bwconncomp(highlight,4);
lenths = [];
    for ic = 1:numel(Ls.PixelIdxList)
        yy = Ls.PixelIdxList{ic};
        if Para.time(max(yy))-Para.time(min(yy)) < 0.1
            highlight(Ls.PixelIdxList{ic}) = 0;
        end 
    end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;


hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.time,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'k*');
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);
title(filename(1:end-4))
xlabel('Time relative to camera change (sec)')
ylabel('Normalised Power (a.u)')
ylim([min(y2plot(:))-0.3*range(y2plot(:)),max(y2plot(:))+0.4*range(y2plot(:))]);
plot([0,0],get(gca,'ylim'),'k--')




%% plot PLV
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/PLV/*.mat']);

if ~filterindex 
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end 
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.time,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'k*');
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.01 (corrected)']);
title(filename(1:end-4))
ylabel('Normalised Power (a.u)')
ylim([min(y2plot(:))-0.3*range(y2plot(:)),max(y2plot(:))+0.4*range(y2plot(:))]);


%% plot coherence

clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/COH/*.mat']);

if ~filterindex 
    return
end

% load in TFR LMEM
load([pathname filename])

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end 
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.freq, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.freq, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.freq,(max(y2plot(:))+0.3*range(y2plot(:)))*highlight,'k*');
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);
title(filename(1:end-4))
xlabel('Frequency (Hz)')
ylabel('Coherence (a.u)')
ylim([min(y2plot(:))-0.3*range(y2plot(:)),max(y2plot(:))+0.4*range(y2plot(:))]);



%% plot multivariate granger causality
 
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/Granger/*.mat']);

if ~filterindex 
    return
end

% load in TFR LMEM
load([pathname filename])
% y2plot = y2raw;

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
% [p_fdr, p_masked] = fdr(pMap, 0.05);
% highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end 
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.freq, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.freq, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.freq,(max(y2plot(:))+0.1*range(y2plot(:)))*highlight,'k*');
plot([0,120],[0 0],'k--')
legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);
title(filename(1:end-4))
ti = find(filename=='_');
ylabel([filename(1:ti-1) ' to Granger Index from ' filename(1:ti-1)])

xlim([1 120])
ylim([-1.5*max(abs(y2plot(:))),1.5*max(abs(y2plot(:)))])
      

      
%% plot PSI
clear;
[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/PSI/*.mat']);

if ~filterindex 
    return
end

% load in TFR LMEM
load([pathname filename])
y2plot = yraw;
tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

% uncorrected
% highlight = double(pMap< 0.05);

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05);
highlight = p_masked;

% remove significant clusters shorter than 100ms
% Ls = bwconncomp(highlight,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         yy = Ls.PixelIdxList{ic};
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             highlight(Ls.PixelIdxList{ic}) = 0;
%         end 
%     end
highlight =double(highlight);
highlight(highlight==0) = nan;

% the significant points could be marked with a red stars

hf = figure;

hold on
hM = shadedErrorBar(Para.time, y2plot(1,:),se2plot(1,:),{'color',[255 106 106]/255},1);
hS = shadedErrorBar(Para.time, y2plot(2,:), se2plot(2,:),{'color',[30 144 255]/255},1);
hsig = plot(Para.time,(max(y2plot(:))+0.1*range(y2plot(:)))*highlight,'k*');
title(filename(1:end-4))
xlim([-0.5 1])
ylim([-1.5*max(abs(y2plot(:))),1.5*max(abs(y2plot(:)))])
ti = find(filename=='_');
xlabel('Time relative to camera change (sec)')
ylabel([filename(1:ti-1) ' to  PSI  from ' filename(1:ti-1)])
plot([0,0],get(gca,'ylim'),'k--')

legend([hM.mainLine,hS.mainLine,hsig], ...
    ['Intact'],['Scrambled'],['P<0.05 (corrected)']);
      
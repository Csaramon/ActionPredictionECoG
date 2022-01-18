
clear;

[filename, pathname, filterindex] = uigetfile(['/Users/qinchaoyi/Desktop/ActionPrediction/Results/TFR/*.mat']);

if ~filterindex
    return
else
    % result data
    load([pathname filename])
end




%% plot time frequency power

if contains(pathname,'TFR')


tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

cdat = tMap;
vdat = Para.freq;
hdat = Para.time;
% vdat = Para.freqvec_amp;
% hdat = Para.freqvec_ph;
clim = [min(cdat(:)) max(cdat(:))];

% uncorrected
% highlight = pMap< 0.05;

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05,'Parametric');
highlight = p_masked;

% the significant voxels could be outlined with a black contour
% plot outline
hf = figure;
% subplot(1,3,3)
h = pcolor(hdat,vdat,cdat);shading interp
[x,y] = meshgrid(hdat, vdat);
x = interp2(x, 2); % change to 4 for round corners
y = interp2(y, 2); % change to 4 for round corners
contourlines = highlight==1;

Ls = bwconncomp(contourlines,4);
lenths = [];
    for ic = 1:numel(Ls.PixelIdxList)
        [xx,yy] = ind2sub(size(contourlines),Ls.PixelIdxList{ic});
        if Para.time(max(yy))-Para.time(min(yy)) < 0.1
            contourlines(Ls.PixelIdxList{ic}) = 0;

        end
        
    end

contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
dx = mean(diff(x(1, :))); % remove for round corners
dy = mean(diff(y(:, 1))); % remove for round corners
hold on
contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[238 92 66]/255,'LineWidth',2);
plot([0,0],[min(Para.freq),max(Para.freq)],'k--')

% plotting parameters
caxis([-4 4])
title([filename(1:end-4) ' (Nsub=' num2str(numel(unique(lmeTBL.Sub))) ' Nelec=' num2str(numel(unique(lmeTBL.Elec))) ')'])
xlabel('Time relative to camera change (sec)')
ylabel('Frequency (Hz)')
ch = colorbar('peer',gca,'EastOutside');
set(ch,'position',[0.92 0.12 0.03 0.8],'ticks',[-4,-3,3,4],'fontsize',12, ...
    'ticklabels',{'-4','-3','3','4'})

set(ch.Label,'string',['S>I    t-value    I>S'],...
    'position',[0.9 0.3],'fontsize',12);

% y2plot=y2plot(:,1:32,:);
% a = repmat(mean(y2plot,3),1,1,size(y2plot,3));
% subplot(1,3,1)
% hM = pcolor(hdat,vdat,squeeze((y2plot(1,:,:)-a(1,:,:))./a(1,:,:)));shading interp
% ca = get(gca,'CLim');
% 
% 
% 
% subplot(1,3,2)
% hS = pcolor(hdat,vdat,squeeze((y2plot(2,:,:)-a(2,:,:))./a(2,:,:)));shading interp
% caxis(ca)

end

%% plot time frequency PAC

if contains(pathname,'PAC')


tMap = tMap';
pMap = pMap';
% tMap(tMap==0)=nan;
% pMap(pMap==0)=nan;
% ynorm = yraw./(repmat(mean(yraw,3),1,1,size(yraw,3)));
ynorm = y2plot;
cdat = tMap;
vdat = Para.freqhigh;
hdat = Para.freqlow;
clim = [min(cdat(:)) max(cdat(:))];

% uncorrected
% highlight = pMap< 0.05;

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05,'nonParametric');
highlight = p_masked;

% the significant voxels could be outlined with a black contour
% plot outline
hf = figure;
subplot(1,3,3)
h = pcolor(hdat,vdat,cdat);shading interp
[x,y] = meshgrid(hdat, vdat);
x = interp2(x, 2); % change to 4 for round corners
y = interp2(y, 2); % change to 4 for round corners
contourlines = highlight==1;

% Ls = bwconncomp(contourlines,4);
% lenths = [];
%     for ic = 1:numel(Ls.PixelIdxList)
%         [xx,yy] = ind2sub(size(contourlines),Ls.PixelIdxList{ic});
%         if Para.time(max(yy))-Para.time(min(yy)) < 0.1
%             contourlines(Ls.PixelIdxList{ic}) = 0;
% 
%         end
%         
%     end

contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
dx = mean(diff(x(1, :))); % remove for round corners
dy = mean(diff(y(:, 1))); % remove for round corners
hold on
contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[238 92 66]/255,'LineWidth',2);

% plotting parameters
% caxis([-4 4])
title([filename(1:end-4) ' (Nsub=' num2str(numel(unique(lmeTBL.Sub))) ' Nelec=' num2str(numel(unique(lmeTBL.Elec))) ')'])
ti = find(filename=='_');
xlabel([filename(1:ti-1) ' Phase (Hz)'])
ylabel([filename(ti:end) ' Amplitude (Hz)'])
xlim([4 30])

subplot(1,3,1)
hM = pcolor(hdat,vdat,squeeze(ynorm(1,:,:))');shading interp
ca = get(gca,'CLim');
xlim([4 30])
subplot(1,3,2)
hS = pcolor(hdat,vdat,squeeze(ynorm(2,:,:))');shading interp
caxis(ca)
xlim([4 30])

end

%% plot time frequency Coherence

if contains(pathname,'COHtf') | contains(pathname,'COHtrl')

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;


iind = 1:1:numel(Para.freq);

pMap = pMap(iind,:);
tMap = tMap(iind,:);
y2plot = y2plot(:,iind,:);


cdat = tMap;
vdat = Para.freq(iind);
hdat = Para.timePT;
% hdat = timePt;
clim = [min(cdat(:)) max(cdat(:))];

% uncorrected
% highlight = pMap< 0.05;

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05,'Parametric');
highlight = p_masked;

% the significant voxels could be outlined with a black contour
% plot outline


ymean = (repmat(mean(mean(y2plot,1),3),size(y2plot,1),1,size(y2plot,3)));
y2plot = (y2plot-ymean)./ymean;
hf = figure;

subplot(1,2,1)
hM = pcolor(hdat,vdat,squeeze(y2plot(1,:,:)));shading interp
ca = get(gca,'CLim');
ylabel(['Frequency (Hz)'])
title('Intact')

[x,y] = meshgrid(hdat, vdat);
x = interp2(x, 2); % change to 4 for round corners
y = interp2(y, 2); % change to 4 for round corners
contourlines = highlight==1;

Ls = bwconncomp(contourlines,4);
lenths = [];
    for ic = 1:numel(Ls.PixelIdxList)
        [xx,yy] = ind2sub(size(contourlines),Ls.PixelIdxList{ic});
        if Para.timePT(max(yy))-Para.timePT(min(yy)) < 0.095
            contourlines(Ls.PixelIdxList{ic}) = 0;

        end
        
    end

contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
dx = mean(diff(x(1, :))); % remove for round corners
dy = mean(diff(y(:, 1))); % remove for round corners
hold on
contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[238 92 66]/255,'LineWidth',2);
caxis(ca)

subplot(1,2,2)
hS = pcolor(hdat,vdat,squeeze(y2plot(2,:,:)));shading interp
caxis(ca)
xlabel(['Time (s)'])
title('Scrambled')

ti = find(filename=='_');
ch = colorbar('peer',gca,'EastOutside');
set(ch,'position',[0.92 0.12 0.03 0.8],'ticks',[],'fontsize',12, ...
    'ticklabels',{})
set(ch.Label,'string',['Low      Coherence      High'],...
    'position',[0.9 0],'fontsize',12);

suptitle({'Coherence'; ...
    [filename(1:ti-1) '↔' filename(ti+1:end-4) ' (Nsub=' num2str(numel(unique(lmeTBL.Sub))) ' Nelec=' num2str(numel(unique(lmeTBL.Elec))) ')']});

end

%% plot time frequency Granger Index

if contains(pathname,'GrangerTF')

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

iind = 1:2:119;
% iind = 1:61;

pMap = pMap(iind,:);
tMap = tMap(iind,:);
y2plot = y2plot(:,iind,:);

cdat = tMap;
vdat = Para.freq(iind);
hdat = Para.timePT;

clim = [min(cdat(:)) max(cdat(:))];

% uncorrected
highlight = pMap< 0.05;

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
% [p_fdr, p_masked] = fdr(pMap, 0.05,'Parametric');
% highlight = p_masked;

% the significant voxels could be outlined with a black contour
% plot outline
hf = figure;

subplot(1,2,1)
hM = pcolor(hdat,vdat,squeeze(y2plot(1,:,:)));shading interp
ca = get(gca,'CLim');
ylabel(['Frequency (Hz)'])
title('Intact')

[x,y] = meshgrid(hdat, vdat);
x = interp2(x, 2); % change to 4 for round corners
y = interp2(y, 2); % change to 4 for round corners
contourlines = highlight==1;

Ls = bwconncomp(contourlines,4);
lenths = [];
    for ic = 1:numel(Ls.PixelIdxList)
        [xx,yy] = ind2sub(size(contourlines),Ls.PixelIdxList{ic});
        if Para.timePT(max(yy))-Para.timePT(min(yy)) < 0.095
            contourlines(Ls.PixelIdxList{ic}) = 0;

        end
        
    end

contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
dx = mean(diff(x(1, :))); % remove for round corners
dy = mean(diff(y(:, 1))); % remove for round corners
hold on
contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[238 92 66]/255,'LineWidth',2);
caxis(ca)

subplot(1,2,2)
hS = pcolor(hdat,vdat,squeeze(y2plot(2,:,:)));shading interp
caxis(ca)
xlabel(['Time (s)'])
title('Scrambled')

ti = find(filename=='_');
ch = colorbar('peer',gca,'EastOutside');
set(ch,'position',[0.92 0.12 0.03 0.8],'ticks',[],'fontsize',12, ...
    'ticklabels',{})
set(ch.Label,'string',['← ' filename(ti+1:end-4) '    From    ' filename(1:ti-1) ' →'],...
    'position',[0.9 0],'fontsize',12);

suptitle({'Granger Index'; ...
    [filename(1:ti-1) '↔' filename(ti+1:end-4) ' (Nsub=' num2str(numel(unique(lmeTBL.Sub))) ' Nelec=' num2str(numel(unique(lmeTBL.Elec))) ')']});

end


%% plot time frequency PSI

if contains(pathname,'PSItf')

tMap(tMap==0)=nan;
pMap(pMap==0)=nan;

iind = 1:1:numel(Para.freq);


pMap = pMap(iind,:);
tMap = tMap(iind,:);
y2plot = y2plot(:,iind,:);

cdat = tMap;
vdat = Para.freq(iind);
hdat = Para.timePT;

clim = [min(cdat(:)) max(cdat(:))];

% uncorrected
% highlight = pMap< 0.05;

% Bofforoni  correction
% highlight = pMap< 0.05/numel(pMap); 

% FDR  correction
[p_fdr, p_masked] = fdr(pMap, 0.05,'Parametric');
highlight = p_masked;

% the significant voxels could be outlined with a black contour
% plot outline
hf = figure;

subplot(1,2,1)
hM = pcolor(hdat,vdat,squeeze(yraw(1,:,:)));shading interp
ylim([5 120])
caxis([-0.0002 0.0002])
ylabel(['Frequency (Hz)'])
title('Intact')

[x,y] = meshgrid(hdat, vdat);
x = interp2(x, 2); % change to 4 for round corners
y = interp2(y, 2); % change to 4 for round corners
contourlines = highlight==1;

Ls = bwconncomp(contourlines,4);
lenths = [];
    for ic = 1:numel(Ls.PixelIdxList)
        [xx,yy] = ind2sub(size(contourlines),Ls.PixelIdxList{ic});
        if Para.timePT(max(yy))-Para.timePT(min(yy)) < 0.095
            contourlines(Ls.PixelIdxList{ic}) = 0;

        end
        
    end

contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
dx = mean(diff(x(1, :))); % remove for round corners
dy = mean(diff(y(:, 1))); % remove for round corners
hold on
contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[238 92 66]/255,'LineWidth',2);


subplot(1,2,2)
hS = pcolor(hdat,vdat,squeeze(yraw(2,:,:)));shading interp
ylim([5 120])
caxis([-0.0002 0.0002])
xlabel(['Time (s)'])
title('Scrambled')

ti = find(filename=='_');
ch = colorbar('peer',gca,'EastOutside');
set(ch,'position',[0.92 0.12 0.03 0.8],'ticks',[],'fontsize',12, ...
    'ticklabels',{})
set(ch.Label,'string',['← ' filename(ti+1:end-4) '    From    ' filename(1:ti-1) ' →'],...
    'position',[0.9 0],'fontsize',12);

suptitle({'Phase Slope Index'; ...
    [filename(1:ti-1) '↔' filename(ti+1:end-4) ' (Nsub=' num2str(numel(unique(lmeTBL.Sub))) ' Nelec=' num2str(numel(unique(lmeTBL.Elec))) ')']});

end

% save the figure to data location
saveas(hf,[pathname filename(1:end-4)])
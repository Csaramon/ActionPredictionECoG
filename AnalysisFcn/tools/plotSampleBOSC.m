tt = 258;  % dongfenglian elec97: 111; 195;421

% liuzhiming elec5: 258

time2plot = [tt tt+6];
fsample = 500;
freqIndex = (eBOSCresults.F >= 3.2 & eBOSCresults.F <= 20);
timeIndex = (eBOSCresults.T >= time2plot(1) & eBOSCresults.T <= time2plot(2));
foi = eBOSCresults.F(freqIndex);
toi = eBOSCresults.T(timeIndex);

episodes_b = eBOSCresults.episodes_b;
% choose episode in foi
ie = 1;
while ie <= size(episodes_b,1)
    if episodes_b{ie,3} >= foi(1) & episodes_b{ie,3} <= foi(10) & ...
            (episodes_b{ie,1}(1,2)./fsample) >= time2plot(1) & (episodes_b{ie,1}(end,2)./fsample) <= time2plot(2)
        ie = ie+1;
    else
        episodes_b(ie,:) = [];
    end
end

BOSCresults.B_b =zscore(BOSCresults.B_b,0,2);
tfr = BOSCresults.B_b(freqIndex,timeIndex);

rawSig = ftDataAll.trial{1}(5,timeIndex);

hraw = subplot(2,1,1);
plot(hraw,toi,rawSig,'k','linewidth',1.5);
hold on
ylabel('Amplitude (μV)')
ylim([-100 150])

htfr = subplot(2,1,2);
pcolor(htfr,toi,foi,tfr);
hold on
xlabel('Elapsed time (sec)')
ylabel('Frequency (Hz)')
% set(gca,'yscale','log')
% set(gca,'ytick',[4 8 16 32])
shading interp
hold on
% choose episode in foi
for ie = 1:size(episodes_b,1)
tt = episodes_b{ie,1};
tin = tt(:,2);
tt(:,2) = tt(:,2)./fsample;
plot(hraw,tt(:,2),ftDataAll.trial{1}(5,tin),'r','linewidth',1.5)

plot(htfr,tt(:,2),eBOSCresults.F(tt(:,1)),'r','linewidth',2)
end

eventdata_b = EventOperation(eBOSCresults.eventdata_b,'eBOSC');
true_latency = eventdata_b.true_latency./1000;
marker2plot = true_latency(true_latency>time2plot(1) & ...
    true_latency<time2plot(2));

for im = 1:size(marker2plot,1)
    plot([marker2plot(im) marker2plot(im)],[foi(1) foi(end)],'k','linewidth',1)
end

latency = eventdata_b.latency./1000;
marker2plot = latency(latency>time2plot(1) & ...
    latency<time2plot(2));

for im = 1:size(marker2plot,1)
    plot([marker2plot(im) marker2plot(im)],[foi(1) foi(end)],'w','linewidth',1)
end
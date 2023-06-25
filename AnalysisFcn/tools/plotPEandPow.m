tp  = 1:90;
elec = 2;
figure;
hPE = plot(tp,allPE(tp),'linewidth',1.5);
hold on
% Y = interp(allpow2sav(tp,elec),2);
tp2 = interp(tp,2);
pow2plot = allpow2sav(tp,elec)+abs(min(allpow2sav(tp,elec)));
hPow = plot(tp,pow2plot./max(pow2plot),'linewidth',1.5);
hSwitch = stairs(tp,ally{1}.*0.5+0.25,'linewidth',1,'color',[0.75 0.75 0.75]);
legend([hPE,hPow,hSwitch],['Prediction Error'],['Normalised Power'],['Perceptual Status'])
xlabel('Time course (second)')

set(gca, 'TickDir','out')
box off
set(gca, 'LineWidth',1.5)
xlim([10 70])

% calculate significant roc electrode

% for iroi = 1:numel(roiLabel)
%     disp(roiLabel{iroi})
% numel(find(allROCBNS{iroi}(:,2)<0.05))./length(allROCBNS{iroi})
% 
% % numel(find(allROCBR{iroi}(:,2)<0.05))./length(allROCBR{iroi})
% 
% % numel(find(allROCRRN{iroi}(:,2)<0.05))./length(allROCRRN{iroi})
% end
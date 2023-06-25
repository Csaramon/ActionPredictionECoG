function [peak]=getpeak(Pow,frequency,ll,ul,threshold)
% this is a peak picking program for mutiple peaks
pows = log10(Pow);
freqs=log10(frequency); % transform to log space
brob=robustfit(freqs(freqs>=ll&freqs<=ul),pows(freqs>=ll&freqs<=ul)); %robust fit
sub=pows-(brob(1)+brob(2)*freqs);
TD=mean(sub)+threshold*std(sub);


[~,peakinds] = findpeaks(sub);
peakvec = zeros(1,length(sub));
peakvec(peakinds) = 1;



peak=frequency(sub>TD&peakvec&freqs>ll&freqs<ul);

% peak=frequency(sub>TD&islocalmax(sub)&freqs>ll&freqs<ul); %return peaks satisfying all criteria



% islocalmax


% figure;
% plot(freqs,pows); hold on;
% 
% plot(freqs,brob(1)+brob(2).*freqs);


end





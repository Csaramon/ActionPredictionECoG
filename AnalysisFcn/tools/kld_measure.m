 function mival = kld_measure(fs,phase_sig, amp_sig,ext_win)

% function mival = mi_measure(phase_sig, amp_sig)
%
% Returns a value for the MI measure calculated between two signals.
% (Functionality to deal with multiple trials will be added soon)
%
% INPUTS:
%
% phase_sig - the instantaneous phase values for a signal which has been
% filtered for a lower, modulating frequency, passed as a column vector
%
% amp_sig - the amplitude values for a signal which has been filtered for a
% higher, modulated frequency, passed as a column vector 
%
% Author: Angela Onslow, May 2010


%% original code(every trial have a pac value )

% num_trials = size(phase_sig, 2);
% 
% for count = 1:num_trials
%     
%     %Create composite signal
%     z = amp_sig(:,count).*exp(1i*phase_sig(:,count));
%     
%     z = abs(sum(z))./sum(amp_sig(:,count));
%     
%     m_raw(count) = z;  %Compute the mean length of composite signal.   weighted PAC in 0 ~ 1
%     
% % %     m_raw(count) = mean(z);  %Compute the mean length of composite signal.
% %        
% % %     mival(count,1) = abs((m_raw(count)));
% 
% 
% 
% end
% % 
% if num_trials > 1
%     mival = nanmean(m_raw);
% end



%% effect of epoch filtering on phase and amplitude

% figure;
% 
% plot(amp_sig);


%% all_in_one MI value(all trials joint into one mat and have only one pac over all trials)
amp_sig = amp_sig(ext_win*fs+1:end- ext_win*fs,:);

phase_sig=phase_sig(ext_win*fs+1:end- ext_win*fs,:);


amp_sig_vec = reshape(amp_sig,[size(amp_sig,1)*size(amp_sig,2) 1]);

ph_sig_vec = reshape(phase_sig,[size(phase_sig,1)*size(phase_sig,2) 1]);


yaa = amp_sig_vec;
ypp = ph_sig_vec+pi; 
        
n_bins = 20;
am_bin = zeros(1,n_bins);
for nb=1:n_bins
    ph_bin = find(ypp>=(nb-1)*pi/(n_bins/2) & ypp<nb*pi/(n_bins/2));
    if ~isempty(ph_bin)
        am_bin(nb) = mean(yaa(ph_bin));
    end;
end;

am_bin = am_bin ./ sum(am_bin);
% bar(am_bin);
infcont = zeros(size(am_bin));
for ic=1:length(am_bin)
    if am_bin(ic) == 0
        infcont(ic) = 0;
    else
        infcont(ic) = am_bin(ic)*log(am_bin(ic));
    end;
end;

mival = 1 + (sum(infcont)/log(n_bins));




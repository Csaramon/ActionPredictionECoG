function [pacmat, pacmatSig] = find_pac (sig_pac, sig_mod, chancmb, ph_freq_vec, amp_freq_vec, time_range, shf_time, pvalue, plt)

% default parameters
if nargin < 2
    sig_mod = sig_pac;
end
if nargin < 3
    chancmb = [1 1];
end
if nargin < 4
    ph_freq_vec = [3 30];
end
if nargin < 5
    amp_freq_vec = [32 200];
end
if nargin < 6
    time_range = [-1 1];
end
if nargin < 7
    shf_time = 0;
end
if nargin < 8
    pvalue = 0.05;
end
if nargin < 9
    plt = 'n';
end

% check datatype
[sig_pac] = ft_checkdata(sig_pac,'datatype','freq');
[sig_mod] = ft_checkdata(sig_mod,'datatype','freq');
if isfield(sig_pac,'fourierspctrm') && isfield(sig_mod,'fourierspctrm')
else error('Input data need fourierspctrm field.');end

% extract frequency points and time points
sig_mod_fp = find(sig_mod.freq >= ph_freq_vec(1) & sig_mod.freq <= ph_freq_vec(2));
sig_pac_fp = find(sig_pac.freq >= amp_freq_vec(1) & sig_pac.freq <= amp_freq_vec(2));
timepoints = find(sig_pac.time >= time_range(1) & sig_pac.time <= time_range(2));

if isempty(sig_mod_fp) || isempty(sig_pac_fp)
    error('No specific frequency found!')
end
if isempty(timepoints)
    error('No specific time point found!')
end

pacmat = zeros(size(chancmb,1),numel(sig_pac_fp),numel(sig_mod_fp));
shf_pacmat = zeros(size(chancmb,1),shf_time,numel(sig_pac_fp),numel(sig_mod_fp));
pacmatSig = zeros(size(chancmb,1),numel(sig_pac_fp),numel(sig_mod_fp));
total_loop_time = size(chancmb,1)*numel(sig_mod_fp)*numel(sig_pac_fp);
loop_count = 0;
strlen = 0;
fprintf(newline)

% calculate PAC
for ichan = 1:size(chancmb,1)
    sig_mod_oi = sig_mod.fourierspctrm(:,chancmb(ichan,1),:,timepoints);
    sig_pac_oi = sig_pac.fourierspctrm(:,chancmb(ichan,2),:,timepoints);
    
    for mod_freq = 1:numel(sig_mod_fp)
        for pac_freq = 1:numel(sig_pac_fp)
            
            mod_sig = squeeze(sig_mod_oi(:,:,sig_mod_fp(mod_freq),:));
            amp_sig_ph = squeeze(sig_pac_oi(:,:,sig_mod_fp(mod_freq),:));
            amp_sig = squeeze(sig_pac_oi(:,:,sig_pac_fp(pac_freq),:));
            rho = circular_linear_corr(mod_sig,amp_sig_ph,amp_sig);
            pacmat(ichan,pac_freq,mod_freq) = abs(rho);
            
            % shuffle process
            for ishf = 1:shf_time
                shf_mod_sig = shuffle_sig(mod_sig);
                rho_shf = circular_linear_corr(shf_mod_sig,amp_sig_ph,amp_sig);
                shf_pacmat(ichan,ishf,pac_freq,mod_freq) = abs(rho_shf);
                
            end
            
            % display total progress
            s = ['Current progress :' num2str(round(100*loop_count/total_loop_time)) '%'];
            strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
            strlen = strlentmp - strlen;
            loop_count = loop_count+1;
            
            
        end
    end
    
end

pacmatSig = ones(size(pacmat));

% statistic test
for icmb = 1:size(chancmb,1)
    
    if shf_time ~= 0
        for i = 1:size(pacmat,2)
            for j = 1:size(pacmat,3)
                [h, p] = sig_test(pacmat(icmb,i,j), squeeze(shf_pacmat(icmb,:,i,j)), pvalue);
                if h == 0
                    pacmatSig(icmb,i,j) = 0;
                end
            end
        end
    end
end

% plot pac result
if strcmp(plt,'y')
    fprintf('Drawing PAC comodulogram... \n')
    for iplot = 1:size(chancmb,1)
        
        onepacmatSig = pacmatSig(iplot,:,:);
        to_plot = squeeze(pacmat(iplot,:,:).*onepacmatSig);
        hf = figure('Name',['chan ' num2str(chancmb(iplot,1)) '-' num2str(chancmb(iplot,2))]);
        ha = axes('parent',hf);
        hp = pcolor(ha,to_plot);
        shading interp
        xlabel(ha,'Modulating Frequency/Hz')
        ylabel(ha,'Modulated Frequency/Hz')
        freq_xtick = get(ha,'xtick');
        freq_ytick = get(ha,'ytick');
        set(ha,'xticklabel',round(10*sig_mod.freq(sig_mod_fp(freq_xtick)))/10)
        set(ha,'yticklabel',round(10*sig_pac.freq(sig_pac_fp(freq_ytick)))/10)
        saveas(hf,['PAC/chan ' num2str(chancmb(iplot,1)) '-' num2str(chancmb(iplot,2))],'fig')
%         close(hf)
    end
end



end


% calculate circular linear correlation coefficient
function [rho] = circular_linear_corr(mod_sig,amp_sig_ph,amp_sig)

if size(mod_sig,1) ~= 1
    mod_sig = reshape(mod_sig',1,size(mod_sig,1).*size(mod_sig,2));
    amp_sig_ph = reshape(amp_sig_ph',1,size(amp_sig_ph,1).*size(amp_sig_ph,2));
    amp_sig = reshape(amp_sig',1,size(amp_sig,1).*size(amp_sig,2));
end

if size(mod_sig,1) ~= 1 || size(amp_sig,1)
    trialnum = min(size(mod_sig,1),size(amp_sig,1));
end

rho_mat = zeros(trialnum,1);
for itrl = 1:trialnum
    
    low_phase1 = angle(mod_sig(itrl,:));
    low_phase2 = angle(amp_sig_ph(itrl,:));
    high_amplitude = abs(amp_sig(itrl,:));
    highphase = angle(amp_sig(itrl,:));
    
    % calculate circular linear correlation coefficient
    Rca = corr(cos(low_phase1)',high_amplitude');
    Rsa = corr(sin(low_phase1)',high_amplitude');
    Rcs = corr(sin(low_phase2)',cos(highphase)');
    rho_mat(itrl) =sqrt((Rca.^2+Rsa.^2-2*Rca*Rsa*Rcs)./(1-Rcs.^2));
end

% Fisher's z-transform
Z_rho = 1/2*log((1+rho_mat)./(1-rho_mat));

% average across trials
rho= nanmean(Z_rho);

end


function shf_sig = shuffle_sig(sig)

% sig dimension = trial*time
% random time point shuffle
for i = 1:size(sig,1)
    cut_point = round(((0.1+0.8*rand(size(sig,1),1))*size(sig,2)));
    shf_sig(i,:) = [sig(i,cut_point(i):end) sig(i,1:cut_point(i)-1)];
    
end

end

function [h,p] = sig_test(pacval, shf_pacmat, alpha)

p = mean (shf_pacmat > pacval);

h = (p <= alpha);
end

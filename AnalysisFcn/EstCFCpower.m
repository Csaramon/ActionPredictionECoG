% Downsample Data
%                                 cfg = [];
%                                 cfg.resamplefs = 500;
%                                 cfg.detrend         = 'no';
%                                 [testData] = ft_resampledata(cfg, testData);


PowThresh = 0.7;
numIter = 1000;
allCutOffFreq = zeros(numIter,120);
strlen = 0;
for i = 1: numIter
    % display permutation progress
    s = ['Permutation progress: ' num2str(i) '/' num2str(numIter)];
    strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
    strlen = strlentmp - strlen;
    
    
    testData.trial{1} = randn(1,501);
    
    % time-frequency decomposition (multi-taper)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.foi          = 1:120;
    cfg.t_ftimwin = 5./cfg.foi;
    cfg.taper      =  'hanning';
    cfg.toi          = 'all';
    %                                 cfg.precision = 'single';
    % cfg.pad='nextpow2';
    ft_info off
    ft_debug off
    ft_notice off
    ft_warning off
    
    freq = ft_freqanalysis(cfg,testData);
    
    cutOffFreq = zeros(size(freq.freq));
    lpFreq = 4;
    for hFreq = find(freq.freq >= 50 & freq.freq <= 100)
        HFdata = squeeze(freq.powspctrm(:,hFreq,:));
        HFdata(isnan(HFdata)) = [];
        L = numel(HFdata);
        Y = fft(HFdata);
        P2 = abs(Y/L);
        P1HF = P2(1:L/2+1);
        
        % find the cut off frequency for low pass filtering
        while cutOffFreq(hFreq) == 0
            lpHFdata = lowpass(HFdata,lpFreq,testData.fsample,'Steepness',0.99);
            Y = fft(lpHFdata);
            P2 = abs(Y/L);
            P1lpHF = P2(1:L/2+1);
            
            if mean(P1lpHF) < PowThresh*mean(P1HF)
                lpFreq = lpFreq + 0.1;
            else
                cutOffFreq(hFreq) = lpFreq-0.1;
            end
        end
        
        
    end
    allCutOffFreq(i,:) = cutOffFreq;
    
end

%% test fft with plot

% Fs = 500;            % Sampling frequency
% L = numel(HFdata);
% Y = fft(HFdata);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% f = Fs*(0:(L/2))/L;
% plot(f,P1)
% mean(P1)
% hold on
% Y2 = fft(lpHFdata);
% P2 = abs(Y2/L);
% P1 = P2(1:L/2+1);
% f = Fs*(0:(L/2))/L;
% plot(f,P1)
% mean(P1)

% Downsample Data
%                                 cfg = [];
%                                 cfg.resamplefs = 500;
%                                 cfg.detrend         = 'no';
%                                 [testData] = ft_resampledata(cfg, testData);


PowThresh = 0.7;
numIter = 1;
allCutOffFreq = zeros(numIter,120);
strlen = 0;
for i = 1: numIter
    
    testData.trial{1} = randn(1,1001);
    
    % time-frequency decomposition (multi-taper)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.foi          = 1:120;
    cfg.t_ftimwin = 3./cfg.foi;
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
    warning('off','all')
    for hFreq = find(freq.freq >= 30 & freq.freq <= 120)
        lpFreq = 4;
        HFdata = squeeze(freq.powspctrm(:,hFreq,:));
        HFdata(isnan(HFdata)) = [];
        L = numel(HFdata);
        f = testData.fsample*(0:(L/2))/L;
        Y = fft(HFdata);
        P2 = abs(Y/L).^2;
        P1HF = P2(1:L/2+1);
        f(1) = [];
        P1HF(1) = [];
        
        % find the cut off frequency for low pass filtering
        while cutOffFreq(hFreq) == 0
            
                % display permutation progress
    s = ['Permutation progress: ' num2str(i) '/' num2str(numIter) ';Power freq: ' num2str(freq.freq(hFreq)) '; Lowpass freq: ' num2str(lpFreq) 'Hz'];
%     s = ['Permutation progress: ' num2str(i) '/' num2str(numIter)];
    
    strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], s);
    strlen = strlentmp - strlen;
%     
%             lpHFdata = lowpass(HFdata,lpFreq,testData.fsample,'Steepness',0.9);
%             Y = fft(lpHFdata);
%             P2 = abs(Y/L).^2;
%             P1lpHF = P2(1:L/2+1);
%             P1lpHF(1) = [];
            
%             if mean(P1lpHF) < PowThresh*mean(P1HF)
% 
%                 lpFreq = lpFreq + 0.5;
%             else
%                 cutOffFreq(hFreq) = lpFreq;
%             end

            if PowThresh*sum(P1HF) > sum(P1HF(1:lpFreq))
                lpFreq = lpFreq + 1;
            else
                cutOffFreq(hFreq) = f(lpFreq-1);
            end


        end   
    end
    allCutOffFreq(i,:) = cutOffFreq;
    
end

plot(mean(allCutOffFreq),freq.freq)
%% test fft with plot

% Fs = 1000;            % Sampling frequency
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
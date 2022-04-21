% testData.trial{1} = randn(1,1001);
t = testData.time{1};
figure;plot(t,testData.trial{1})

% time-frequency decomposition (multi-taper)
    cfg              = [];
    cfg.output       = 'fourier';
    cfg.method       = 'mtmconvol';
    cfg.foi          = 1:120;
    cfg.t_ftimwin = 2./cfg.foi;
    cfg.taper      =  'hanning';
    cfg.toi          = 'all';
    %                                 cfg.precision = 'single';
    % cfg.pad='nextpow2';
    ft_info off
    ft_debug off
    ft_notice off
    ft_warning off
    
    freq = ft_freqanalysis(cfg,testData);
    
    LFdata = squeeze(freq.fourierspctrm(:,:,20,:));
    HFdata = abs(squeeze(freq.fourierspctrm(:,:,80,:))).^2;

    
    figure();
    subplot(2,1,1);
        plot(t,real(LFdata))
        ylabel('Amplitude')
        legend('20Hz')
        subplot(2,1,2);
        
        
    plot(t,HFdata)
    legend('80Hz')
    ylabel('Power')
    xlabel('Time (sec)')
    title('2 cycles')
    
    
        HFdata(isnan(HFdata)) = [];
    L = numel(HFdata);
    f = testData.fsample*(0:(L/2))/L;
    Y = fft(HFdata);
    P2 = abs(Y/L).^2;
    P1HF = P2(1:L/2+1);
    f(1) = [];
    P1HF(1) = [];
    figure();plot(f,P1HF)
    xlim([0 120])
    ylabel('Power')
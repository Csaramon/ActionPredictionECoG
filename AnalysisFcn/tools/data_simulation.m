cfg = [];
cfg.method     = 'phalow_amphigh';
cfg.fsample    = 2000;
cfg.trllen     = 10;
cfg.numtrl     = 20;
cfg.output     = 'all';
% amplitude modulation
cfg.s1.freq    = 6;
cfg.s1.phase   = 'random'; %phase differs over trials
cfg.s1.ampl    = 1;
% frequency that becomes modulated
cfg.s2.freq    = 90;
cfg.s2.phase   = 'random'; %phase differs over trials
cfg.s2.ampl    = 0.2; 
% DC shift of S1
cfg.s3.freq    = 0;
cfg.s3.phase   = 0;
cfg.s3.ampl    = 1; %determines amount of modulation, should be at least s1.ampl
% noise
cfg.noise.ampl = 0.1;

data = ft_freqsimulation(cfg);
figure;
sel = 1:cfg.fsample;
subplot(3,3,1); plot(data.trial{1}(1,sel)); title(data.label{1})
subplot(3,3,2); plot(data.trial{1}(2,sel)); title(data.label{2})
subplot(3,3,3); plot(data.trial{1}(3,sel)); title(data.label{3})
subplot(3,3,4); plot(data.trial{1}(4,sel)); title(data.label{4})
subplot(3,3,5); plot(data.trial{1}(5,sel)); title(data.label{5})
print -dpng phalow_amphigh_fig1.png

% mix low frequency wave and high frequency PAC wave
% for itrl = 1:numel(data.trial)
%     data.trial{itrl}(1,:) = data.trial{itrl}(1,:)+data.trial{itrl}(2,:);
% end
% show powerspectrum simulated data
% cfg = [];
% cfg.method    = 'mtmfft';
% cfg.channel   = 'mix';
% cfg.output    = 'pow';
% cfg.taper     = 'hanning';
% cfg.foilim    = [2 80];
% cfg              = [];
% cfg.output       = 'fourier';
% cfg.method       = 'hilbert';
% cfg.foi          = logspace(log10(4),log10(180),41);
% cfg.width        =  0.1*cfg.foi; % fractional bandpass width
%                 cfg.toi          =  -1:0.005:1;
% cfg.toi = 'all';
% cfg.keeptrials   = 'yes';
% cfg.filttype = 'firws';
% cfg.filtdir = nan;
% cfg.filtorder = nan;
% cfg.precision = 'single'; % reduce data size
%         cfg.pad = 3;
% ft_spect         = ft_freqanalysis(cfg,data);

% fft_data = ft_freqanalysis(cfg,data);
% figure; ft_singleplotER([],fft_data);
% print -dpng phalow_amphigh_fig2.png
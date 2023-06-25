figure
duration = 3;
interval = 0;
freq = 1000;
fsample = 44100;
repetition = 6;

x = dataIn;
N = length(x);
n = 0:N-1;
t = n/fsample;
% x = sin(2*pi*freq*t);


% xo = dataOut;
% xo2 = dataOut;
% Filtering
% [b,a] = BUTTER(4,2000./fix(fsample./2),'low');
% x = filter(b,a,x);
% x = medfilt1(x,16);

y = fft(x,N);
mag = abs(y);
f = n*fsample/N;
plot(f,mag);

% plot(x)

xlim([0 200])  
xlabel('Frequency/Hz')
ylabel('Amplitude')


dB_Gain = -36;  %-32 equals to 70dB at 1m (amplifier at max-5)


sf = 10.^(dB_Gain./20);
pause(1)
sound(sf.*x,fsample)

pause(1)
sound(sf.*dataOut1,fsample)

pause(1)
sound(sf.*dataOut2,fsample)

pause(1)
sound(sf.*dataOut3,fsample)

pause(1)
sound(sf.*dataOut4,fsample)

pause(1)
sound(sf.*dataOut5,fsample)

bsFilt = designfilt('bandstopfir', ...
    'StopbandFrequency1',4500, ...
    'StopbandFrequency2',5000, ...
    'StopbandAttenuation',3)

fvtool(Hd)
dataIn = randn(22050,1);
dataOut1 = filter(bs1,dataIn);   
dataOut2 = filter(bs2,dataIn); 
dataOut3 = filter(bs3,dataIn); 
dataOut4 = filter(bandstop4,dataIn); 
dataOut5 = filter(bandstop5,dataIn); 


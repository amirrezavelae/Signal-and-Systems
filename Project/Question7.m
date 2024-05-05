clear all;
clc;

Data = importdata("400101967.mat");
Fs = Data.fs;
sigNoise = Data.corrupted_signal;


%% removing signal using LMS algorithm
clear sound;
clc;

% define error noise estimation and filter coeffs

freq = 50;

t = 0 : 1/Fs : (length(sigNoise)-1)/Fs;
refNoise = sin(2*pi*freq * t) + cos(2*pi*freq * t);

step = 0.001;
weiLen = 50;

WInit = rand(1,weiLen);

for k = 1:10
    for i = weiLen:length(t)
        noiseEst(i) = WInit * refNoise(i:-1:i-50+1)';
        sigEst(i) = sigNoise(i) - noiseEst(i);
        WInit = WInit + step * sigEst(i) * refNoise(i:-1:i-50+1);
    end
end

subplot(3,1,1);
plot(t,sigNoise);
title('Noisy Signal','Interpreter','latex');
ylim([-1 1]);

subplot(3,1,2);
plot(t,sigEst);
title('Signal Estimaion','Interpreter','latex');
ylim([-1 1]);

subplot(3,1,3);
plot(t,sigNoise - sigEst);
title('Noise Estimation','Interpreter','latex');
ylim([-1 1]);

% sound(sigEst, Fs);

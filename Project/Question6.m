% Ehsan Merrikhi - Amirreza Vellai
% https://youtu.be/3dLHlxtLy_w
% thanks to our foreign friend for helping us in thie project

%% Question 1
clear all;
% recording voice
% -------------------------------------------------

% setting properties of voice
duration = 5;
Fs = 10^4;
nBits = 24;

% creating record object
s = audiorecorder(Fs, nBits, 1);
recordblocking(s, duration);
% -------------------------------------------------

audiowrite("audio.wav",s.getaudiodata, Fs);

%% adding 50Hz noise
clc;
clear all;


[sig,Fs] = audioread("audio.wav");
% -------------------------------------------------
strPhs = rand * 2*pi;
amp = 0.5;
freq = 50;

t = 0 : 1/Fs : 5-1/Fs;
noise = amp * sin(t * 2*pi * freq + strPhs);

sigNoise = sig + noise';
% -------------------------------------------------

% voice = audioplayer(sigNoise, Fs);
% play(voice);

%% removing signal using LMS algorithm
clc;
clear("noise");

freq = 50;


t = 0 : 1/Fs : 5-1/Fs;
refNoise = sin(2*pi*freq * t) + cos(2*pi*freq * t);

step = 0.0005;
weiLen = 50;

% initiating weight function
WInit = rand(1,weiLen);

% applying alorithm
for k = 1:10
    for i = weiLen:length(t)
        noiseEst(i) = WInit * refNoise(i:-1:i-50+1)';
        sigEst(i) = sigNoise(i) - noiseEst(i);
        WInit = WInit + step * sigEst(i) * refNoise(i:-1:i-50+1);
    end
end

subplot(4,1,1);
plot(t,sig);
title('Original Signal','Interpreter','latex');
ylim([-1 1]);

subplot(4,1,2);
plot(t,sigNoise);
title('Noisy Signal','Interpreter','latex');
ylim([-1 1]);

subplot(4,1,3);
plot(t,sigEst);
title('Signal Estimation','Interpreter','latex');
ylim([-1 1]);

subplot(4,1,4);
plot(t,sig - sigEst');
title('Noise Estimation','Interpreter','latex');
ylim([-1 1]);

t = timer;
t.StartDelay = 5;
t.TimerFcn = @(~,~) sound(3*sigNoise, Fs);

sound(3*sigEst, Fs);
start(t);

audiowrite('noisy signal.wav',sigNoise,Fs);
audiowrite('recoverd signal.wav',sigEst,Fs);


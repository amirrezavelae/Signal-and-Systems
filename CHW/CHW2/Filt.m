function Hd = Filt
%FILT Returns a discrete-time filter object.

% Butterworth Bandpass filter

% All frequency values are in Hz.
Fs = 48000;  % Sampling Frequency

N   = 10;
Fc1 = 50;
Fc2 = 1200;

h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

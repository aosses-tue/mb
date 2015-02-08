% Received on 20150206 from AK (e-mail)

% My description:
%       inputfile - 't54.wav'
%
%       inputfile is low-pass filtered. After being filtered, the resulting 
%       sound is played.

clear
clf
inputfile = 't54.wav';
[s1, sampleFrequency] = wavread(inputfile);
t = [1:length(s1)]/sampleFrequency;
cutoffFrequency = 200;
% calculate coefficients of lowpass filter
ord = 4;
[b, a ] = butter(ord, 2*cutoffFrequency/sampleFrequency);
s2 = s1;
fact = 256;
% s2 = round(fact*s1)/fact;
s2 = ord*filter(b, a, s1);
sound(s2, sampleFrequency)

subplot(2,1,1)
plot(t, s1)
axis([0.1*t(length(t)) 0.15*t(length(t)) -1 1])
title(sprintf('%s',inputfile))

subplot(2,1,2)
plot(t, s2)
axis([0.1*t(length(t)) 0.15*t(length(t)) -1 1])
title(sprintf('+ LPF, fc=%.0f [Hz]',cutoffFrequency))

disp('')
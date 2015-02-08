% Received on 20150206 from AK (e-mail)
%
% My description:
%       inputfile - 't54.wav'
%       outputfile - 't54_8.wav'
%
% Possible input parameters: 
%       inputfile
%       nrOfBits
% 
% Num. of bits reduction: number of amplitude points = 2^(bits)
%       if 2 bits:      amplitudes are [-1 -0.5 0 0.5 1]

clear all
close all
inputfile = 't54.wav';
[s1, sampleFrequency] = wavread(inputfile);
t = (1:length(s1))/sampleFrequency;
nrOfBits = 2;
s = round(2^(nrOfBits-1)*s1);
s2 = s/(2^(nrOfBits-1));
sound(s2, sampleFrequency)

outputfile = sprintf('t54_%.0f.wav',nrOfBits);
wavwrite(s2, sampleFrequency, 16, outputfile); % wavwrite(y,Fs,N,filename) % N-num of bits

subplot(2,1,1)
plot(t, s1)
axis([0.1*t(length(t)) 0.15*t(length(t)) -1 1])
subplot(2,1,2)
plot(t, s2)
axis([0.1*t(length(t)) 0.15*t(length(t)) -1 1])

disp('')
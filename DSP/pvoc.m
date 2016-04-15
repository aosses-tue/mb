function y = pvoc(x, r, n)
% function y = pvoc(x, r, n)
%
% 1. Description:
%       y = pvoc(x, r, n)  Time-scale a signal to r times faster with phase 
%           vocoder x is an input sound. n is the FFT size, defaults to 1024.  
%           Calculate the 75%-overlapped STFT (hop size of 25%), squeeze it 
%           by a factor of r, inverse spegram.
%       r - times faster
%       n - window length [samples]
% 
% 2. Stand-alone example:
%       [insig fs] = Wavread(filename);
%       timesfaster = 0.75; % sampling frequency is implicit in 'timesfaster'
%       N       = 4096; % N = 1024 was giving some phase problems
%       outsig  = pvoc(insig,timesfaster,N);
%       sound(insig,fs);
%       pause(1.2*length(insig)/fs)
% 
%       sound(outsig,fs);
% 
% 2000-12-05, 2002-02-13 dpwe@ee.columbia.edu.  Uses pvsample, stft, istft
% $Header: /home/empire6/dpwe/public_html/resources/matlab/pvoc/RCS/pvoc.m,v 1.3 2011/02/08 21:08:39 dpwe Exp $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
  n = 1024;
end

% With hann windowing on both input and output, 
% we need 25% window overlap for smooth reconstruction
hop = n/4;
% Effect of hanns at both ends is a cumulated cos^2 window (for
% r = 1 anyway); need to scale magnitudes by 2/3 for
% identity input/output
%scf = 2/3;
% 2011-02-07: this factor is now included in istft.m
scf = 1.0;

%% Calculate the basic STFT, magnitude scaled
X = scf * pvoc_stft(x', n, n, hop); % n-points FFT, n - window length, hop - hop size

%% Calculate the new timebase samples
[rows, cols] = size(X);
t = 0:r:(cols-2);
% Have to stay two cols off end because (a) counting from zero, and 
% (b) need col n AND col n+1 to interpolate

%% Generate the new spectrogram
X2 = pvsample(X, t, hop);

% Invert to a waveform
y = pvoc_istft(X2, n, n, hop)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function Demo_fluct
% function Demo_fluct
%
% 1. Description:
%   Fluctuation implemented by Chalupper and provided to PsySound team. These
%   functions were renamed adding the prefix 'ch_' to all the functions.
% 
%   original location: "..\tb_Psysound32\AudioAnalysers\@LoudnessCF\private\"
%   original name: demo_fluct.m
%
% To program:
%       testnoise
%       amp - used by am
%       tonheit(should be conversion to MEL scale
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created in    : 2000
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Last update on: 07/08/2014 
% Last use on   : 11/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amsigs = [];

try
    sig = ch_testnoise(70,20,20000,1); % creates a tone
catch
    sig = Create_test_noise(70,20,20000,1);
end

% freqfft(sig,44100,4096)
% sig2 = Create_test_noise(70,800,1200,1);

% [sig fs] = wavread('D:\MATLAB_git\tb_Psysound32\Documentation\Examples\RoughnessTest\sound.wav');

% amsig = sig;
amsig=ch_am(sig,2,40,'d',44100,-pi/2,1); %white noise, 1kHz, d=40dB, fmod=2Hz,peak level=80dB, introduces an amplitude modulation

%calculate loudness time pattern
% out=ch_dlm(amsig);
[N, main_N, spec_N]=ch_dlm(amsig);

%calculate loudness fluctuation
lf=ch_fluct(main_N)

amsigs = [amsigs amsig];

%  sig = testnoise(70,20,20000,1); % we do not have tonheit, invtonheit
%[sig fs] = wavread('D:\MATLAB_git\tb_Psysound32\Documentation\Examples\RoughnessTest\sound.wav');

amsig=ch_am(sig,32,40,'d',44100,-pi/2,1); %white noise, 1kHz, d=40dB, fmod=32Hz,peak level=80dB

%calculate loudness time pattern
[N, main_N, spec_N]=ch_dlm(amsig);

%calculate loudness fluctuation
lf=ch_fluct(main_N)

amsigs = [amsigs amsig];

% sig = testnoise(70,20,20000,1);
amsig=ch_am(sig,4,4,'d',44100,-pi/2,1); %white noise, 1kHz, d=40dB, fmod=2Hz,peak level=80dB

%calculate loudness time pattern
[N, main_N, spec_N]=ch_dlm(amsig);

%calculate loudness fluctuation
lf=ch_fluct(main_N)

amsigs = [amsigs amsig];

amsigs = transpose(amsigs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(3,1,1)
plot(To_dB(abs(amsigs(1,:))))

ha = gca;

subplot(3,1,2)
plot(To_dB(abs(amsigs(2,:))))
ha(end+1) = gca;

subplot(3,1,3)
plot(To_dB(abs(amsigs(3,:))))
ha(end+1) = gca;

linkaxes(ha, 'xy');
ylim([-50-20 -20])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
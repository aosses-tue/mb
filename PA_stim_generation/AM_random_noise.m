function y = AM_random_noise(Finf,Fsup,SPL,dur,Fs,Fmod,Mdept,dBFS)
% function y = AM_random_noise(Finf,Fsup,SPL,dur,Fs,Fmod,Mdept,dBFS)
% 
% 1. Description:
%       Creates one frame of AM-'running' noise at Fs, N
%       The input parameters are cutoff frequencies Finf and Fsup.
% 
%       See also: AM_random_noise_BW.m
% 
% 2.1 Example:
%       Finf = 750;
%       Fsup = 1250;
%       Fmod = 4;
%       dur = 4;
%       Mdept = 1;
%       SPL = 70;
%       Fs = 44100;
%       y = AM_random_noise(Finf,Fsup,SPL,dur,Fs,Fmod,Mdept);
% 
%       % If you want to store the output (Wav file):
%       AM_random_noise(Finf,Fsup,SPL,dur,Fs,Fmod,Mdept);
% 
% 2.2 Example:
%       Fmod = 4;
%       dur = 4;
%       Fs = 44100;
%       Finf = 0;
%       Fsup = Fs/2;
%       Mdept = 1;
%       SPL = 70;
%       y = AM_random_noise(Finf,Fsup,SPL,dur,Fs,Fmod,Mdept);
%       sound(y,Fs); 
% 
% Programmed/adapted by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Original file name: AMnoiseRandom.m (by Dik Hermes)
% Created on    : 16/03/2015
% Last update on: 16/03/2015 % Update this date manually
% Last use on   : 01/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 8
    dBFS = 100; % 0 dBFS <--> 100 dB SPL 
end

if nargin < 7
    Mdept = 0; % no modulation as default
end

if nargin < 6
    Fmod  = 4;
end

if nargin < 5
    Fs  = 40960;
end

if nargin < 4
    dur = 200e-3; % 200 ms, leading to N = 8192 at Fs = 40960 Hz
end

if nargin < 3
    SPL = 70;
end

if nargin < 2
    Fsup = 0;
end

if nargin < 1
    Finf = 0;
end

BW      = Fsup-Finf;
Fc      = (Finf+Fsup)/2;

N       = round(Fs*dur);
Sig     = rand(N,1)-0.5; % non-calibrated white noise
dF		= Fs/N;
wstep	= 2*pi/Fs;

for q=1:N
   Sig(q) =	Sig(q)*(1+(Mdept*sin(wstep*Fmod*q)));
end
Sig     = fft(Sig);
FcLoc	= round(Fc/dF);             % bin-number of Fc (location)
Finf    = FcLoc-round(BW/(2*dF)+1); % bin-number of Finf = Fc - BW/2
Fsup    = FcLoc+round(BW/(2*dF)+1); % bin-number of Fsup = Fc + BW/2
Bpass	= max([Finf 1]):min([Fsup N]);
BPmul	= zeros(N,1);BPmul(Bpass) = 1;

Sig     = real(ifft(Sig.*BPmul)); % Inverse FFT of band-passed filter signal
Amp     = From_dB(-dBFS)*( From_dB(SPL-3)/mean(rms(Sig)) );
Sig     = Amp*Sig;

y		= Sig;

if nargout == 0
    
    filename = [Get_TUe_paths('outputs') sprintf('randomnoise-Fc-%.0f_BW-%.0f_Fmod-%.0f_Mdept-%.0f_SPL-%.0f',Fc,BW,Fmod,Mdept,SPL)];
    Wavwrite(y,Fs,filename);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
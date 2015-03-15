function y = AMnoise(Fc,BW,Fmod,Mdept,SPL)
% function y = AMnoise(Fc,BW,Fmod,Mdept,SPL)
%
% 1. Description:
%       Creates one frame of AM-'frozen' noise at Fs=40960, N=8192
% 
% 2. Example:
%       Fmod = 4;
%       Mdept = 1;
%       BW = 500;
%       Fc = 1000;
%       SPL = 70;
%       y = AMnoise(Fc,BW,Fmod,Mdept,SPL);
% 
% Created by: Dik Hermes
% Received in: August 2014
% Additional comments by: Alejandro Osses
% Last edit on: -
% Last used on: 15/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs      = 40960;
N       = 8192; % 200 ms at Fs = 40960 Hz
Sig     = zeros(N,1);
dF		= Fs/N;
wstep	= 2*pi/Fs;
load white.mat; % it loads variable whitenoise
for q=1:N
   Sig(q,1)	= whitenoise(q,1)*(1+(Mdept*sin(wstep*Fmod*q)));
end
Sig     = fft(Sig);
FcLoc	= round(Fc/dF);             % bin-number of Fc (location)
Finf    = FcLoc-round(BW/(2*dF)+1); % bin-number of Finf = Fc - BW/2
Fsup    = FcLoc+round(BW/(2*dF)+1); % bin-number of Fsup = Fc + BW/2
Bpass	= max([Finf 1]):min([Fsup N]);
BPmul	= zeros(N,1);BPmul(Bpass) = 1;

Sig     = real(ifft(Sig.*BPmul)); % Inverse FFT of band-passed filter signal
Sig     = Sig*db2amp(SPL-83)/mean(rms(Sig));

y		= Sig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
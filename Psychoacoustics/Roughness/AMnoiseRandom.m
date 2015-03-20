function y = AMnoiseRandom(Fc,BW,Fmod,Mdept,SPL)
% function y = AMnoiseRandom(Fc,BW,Fmod,Mdept,SPL)
% 
% 1. Description:
%       Creates one frame of AM-'running' noise at Fs=40960, N=8192
% 
% 2. Example:     
%       Fmod = 4;
%       Mdept = 1;
%       BW  = 500;
%       Fc  = 1000;
%       SPL = 70;
%       y = AMnoiseRandom(Fc,BW,Fmod,Mdept,SPL);
% 
% Created by: Dik Hermes
% Received in: August 2014
% Additional comments by: Alejandro Osses
% Last edit on: -
% Last used on: 15/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs      = 40960;
N       = 8192; % 200 ms at Fs = 40960 Hz
Sig     = rand(N,1)-0.5;
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
Sig     = Sig*db2amp(SPL-83)/mean(rms(Sig));

y		= Sig;

if nargout == 0
    
    filename = [Get_TUe_paths('outputs') sprintf('randomnoise-Fc-%.0f_BW-%.0f_Fmod-%.0f_Mdept-%.0f_SPL-%.0f',Fc,BW,Fmod,Mdept,SPL)];
    Wavwrite(y,Fs,filename);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
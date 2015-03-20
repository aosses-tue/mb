function y=createFM(Fc,Fmod,dF,SPL,Fs,N)
% function y=createFM(Fc,Fmod,dF,SPL,Fs,N)	
% 
% 1. Description:
%       Creates an Frequency-Modulated sine tone with central frequency Fc,
%       modulation frequency Fmod, frequency deviation dF, Sound Pressure 
%       Level SPL, sampling frequency Fs and length N (in samples).
% 
% 2. Example:
%       Fc  = 1000;
%       Fmod = 4;
%       Fs  = 40960;
%       N   = 8192;
%       SPL = 70;
%       dF  = 700;
%       y=createFM(Fc,Fmod,dF,SPL,Fs,N);
% 
% 3. Additional information:
%       Tested cross-platform: Yes
% 
% Created by: Dik Hermes
% Received in: August 2014
% Additional comments by: Alejandro Osses
% Last edit on: -
% Last used on: 15/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
	Fs	=	40960;
	N	=	8192;
end

Amp     = db2amp(SPL-80);	% -20 dBFS	<--> 60 dB SPL
wstep	= 2*pi/Fs;
Sig     = zeros(N,1);
dw		= 0;

for q	=	1:1:N;
    dw      = dw + wstep*( Fc + dF*sin(wstep*Fmod*q) );
	Sig(q)	= Amp*sin(dw);
end

y	=	Sig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
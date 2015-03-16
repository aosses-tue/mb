function y = createAM(Fc,Fmod,Mdept,SPL,Fs,N)	
% function y = createAM(Fc,Fmod,Mdept,SPL,Fs,N)	
% 
% 1. Description:
%       Creates an Amplitude-Modulated sine tone with central frequency Fc,
%       modulation frequency Fmod, modulation index Mdept, Sound Pressure 
%       Level SPL, sampling frequency Fs and length N (in samples).
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

wstep	=	2*pi/Fs;
Sig     =	zeros(N,1);

for q	=	1:1:N;
	Sig(q)	=	(1+(Mdept*sin(wstep*Fmod*q)))*sin(wstep*Fc*q);
end

Sig	=	Sig*db2amp(SPL-83)/mean(rms(Sig));
y	=	Sig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
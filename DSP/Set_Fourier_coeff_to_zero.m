function inoutsig = Set_Fourier_coeff_to_zero(inoutsig,fs,N,fmin,fmax)
% function inoutsig = Set_Fourier_coeff_to_zero(inoutsig,fs,N,fmin,fmax)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/08/2015
% Last update on: 06/08/2015 
% Last use on   : 06/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    N       = length(inoutsig);
end

fc      = (fmin + fmax)/2;
BW      = fmax - fmin;
dF		= fs/N;
% wstep	= 2*pi/fs;

inoutsig   = fft(inoutsig);
FcLoc	= round(fc/dF);             % bin-number of Fc (location)
Finf    = FcLoc-round(BW/(2*dF)+1); % bin-number of Finf = Fc - BW/2
Fsup    = FcLoc+round(BW/(2*dF)+1); % bin-number of Fsup = Fc + BW/2
Bpass	= max([Finf 1]):min([Fsup N]);
BPmul	= zeros(N,1);BPmul(Bpass) = 1;
inoutsig   = real(ifft(inoutsig.*BPmul)); % Inverse FFT of band-passed filter signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

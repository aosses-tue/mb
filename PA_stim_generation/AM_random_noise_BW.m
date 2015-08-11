function [SigAM filename Env SigBBN] = AM_random_noise_BW(fc,BW,SPL,dur,fs,Fmod,Mdept,dBFS)
% function [SigAM filename Env SigBBN] = AM_random_noise_BW(fc,BW,SPL,dur,fs,Fmod,Mdept,dBFS)
% 
% 1. Description:
%       Creates one frame of AM-'running' noise at fs, N.
%       The input parameters are centre frequency fc and bandwidth BW
%           
%       The Envelope is applied before the band-pass limitation.
% 
%       See also: AM_random_noise.m
% 
% 2.1 Example:
%       Finf = 750;
%       Fsup = 1250;
%       BW = Fsup - Finf;
%       Fc = (Finf + Fsup)/2;
%       Fmod = 4;
%       dur = 4;
%       Mdept = 1;
%       SPL = 70;
%       Fs = 44100;
%       [y env] = AM_random_noise_BW(Fc,BW,SPL,dur,Fs,Fmod,Mdept);
% 
%       % If you want to store the output (Wav file):
%       AM_random_noise_BW(Fc,BW,SPL,dur,Fs,Fmod,Mdept);
% 
% 2.2 Example:
%       Fmod = 4;
%       dur = 4;
%       Fs = 44100;
%       Finf = 0;
%       Fsup = Fs/2;
%       BW = Fsup - Finf;
%       Fc = (Finf + Fsup)/2;
%       Mdept = 1;
%       SPL = 70;
%       y = AM_random_noise_BW(Fc,BW,SPL,dur,Fs,Fmod,Mdept);
%       sound(y,Fs); 
% 
% Programmed/adapted by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Original file name: AMnoiseRandom.m (by Dik Hermes)
% Created on    : 16/03/2015
% Last update on: 16/03/2015 % Update this date manually
% Last use on   : 16/03/2015 % Update this date manually
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
    fs  = 40960;
end

if nargin < 4
    dur = 200e-3; % 200 ms, leading to N = 8192 at Fs = 40960 Hz
end

if nargin < 3
    SPL = 70;
end

if nargin < 2
    BW  = 100;
end

if nargin < 1
    fc = 1000;
end

Fsup    = fc + BW/2;
Finf    = fc - BW/2;

N       = round(fs*dur);
Sig     = rand(N,1)-0.5; % non-calibrated white noise, with amplitudes from -0.5 to 0.5
SigAM   = zeros(N,1); % Initialisation
Env     = zeros(N,1); % Initialisation
dF		= fs/N;
wstep	= 2*pi/fs;

for i=1:N
    Env(i) = (1+(Mdept*sin(wstep*Fmod*i)));
    SigAM(i) =	Sig(i)*Env(i);
end

SigBBN  = SigAM;
SigAM   = fft(SigAM);
FcLoc	= round(fc/dF);             % bin-number of Fc (location)
Finf    = FcLoc-round(BW/(2*dF)+1); % bin-number of Finf = Fc - BW/2
Fsup    = FcLoc+round(BW/(2*dF)+1); % bin-number of Fsup = Fc + BW/2
Bpass	= max([Finf 1]):min([Fsup N]);
BPmul	= zeros(N,1);BPmul(Bpass) = 1;
SigAM   = real(ifft(SigAM.*BPmul)); % Inverse FFT of band-passed filter signal

SigAM   = setdbspl(SigAM ,SPL,dBFS);
SigBBN  = setdbspl(SigBBN,SPL,dBFS);
% filename will be stored in case only if nargout == 0
filename{1} = sprintf('randomnoise-fc-%.0f_BW-%.0f_fmod-%.0f_Mdept-%.0f_SPL-%.0f',fc,BW,Fmod,100*Mdept,SPL);
filename{2} = sprintf('randomnoise-BBN_SPL-%.0f',SPL);
fullfilename  = [Get_TUe_paths('outputs') filename{1}];
fullfilename2 = [Get_TUe_paths('outputs') filename{2}];

if nargout == 0
    
    Wavwrite(SigAM ,fs,fullfilename);
    Wavwrite(SigBBN,fs,fullfilename2);
    
    figure;
    t = (1:N)/fs;
    plot(t,Env,'r--',t,SigAM,'b');
    legend(sprintf('Envelope that was applied\n before setting coeff to 0'))
    xlabel('Time [s]')
    ylabel('Amplitude')
    grid on
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
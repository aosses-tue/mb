function [y ydB f, RMS] = envfreqfft(x,K,fs,windowtype,dBFS,cutoffLPF)
% function [y ydB f, RMS] = envfreqfft(x,K,fs,windowtype,dBFS,cutoffLPF)
%
% 1. Description:
%       FFT with plot (if no output is specified)
% 
% 3. Stand-alone example:
%       
%       [x fs] = Wavread(filename);
%       windowtype = 'hanning';
%       dBFS = 100;
%       figure;
%       freqfft2(x,K,fs,windowtype,dBFS);     
% 
% Programmed by Alejandro Osses, TUe
% Created on    : 12/08/2015
% Last update on: 12/08/2015
% Last use on   : 12/08/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    cutoffLPF = 1000;   % 1 = semilogx
                        % 2 = linear scale, normalised scale
                        % 3 = linear scale, 0 to fs/2 (if fs is defined)
                        % 4 = FFT bins
end

if nargin < 5
    dBFS = 100;
end

N = length(x);
f = (1:K)/K*fs/2;

if nargin < 4
    windowtype = 'hanning';
end

f = f(:);

h = Get_window(windowtype,N);

[b,a] = butter(2,cutoffLPF/(fs/2));
x = abs(hilbert(x));
x = filter(b,a,x);

switch windowtype
    case 'rectangular'
        corrWindow_dB = 0;
    otherwise
        corr = length(h)*mean(h);
        corrWindow_dB = log10(corr) + 0.272;
end

y   = fft(h.*x,K*2);
y   = y(1:K,:);

df  = f(2) - f(1);

constant = sqrt(1/(df*K*N));
y = From_dB(dBFS) * y * constant * From_dB(corrWindow_dB);

ydB = 20*log10(abs(y));

RMS = rmsdb_freqdomain(f,y);

if nargout == 0
    
    plot(f,ydB), grid on
    xlabel('Frequency (Hz)')
    xlim([20 fs/2])
    title(sprintf('Avg. level = %.2f dB',RMS))

    ylabel('Magnitude [dB]')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
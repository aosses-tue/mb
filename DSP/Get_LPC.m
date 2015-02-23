function [ffreq out] = Get_LPC(x,fs,N)
% function [ffreq out] = Get_LPC(x,fs,N)
%
% 1. Description:
%       Gets the LPC coefficients for the signal x, (sampling freq. of fs).
%       The number of LPC coefficients is determined using the rule of the
%       thumb Ncoeffs = 2+round(fs/1000);
% 
%       Inputs:
%           x  - input signal
%           fs - sampling frequency of the input signal
%           N  - number of points for FFT plot
%       Outputs:
%           out - struct containing t, f, ncoeff, h_dB
% 
% 2. Stand-alone example:
%       filename= 'D:\Output\tmp-Audio\meas-ac-mode-5.wav';
%       [x fs]  = Wavread(filename);
%       N = 4096; 
%       Get_LPC(x,fs,N);
%
%       filename= '~/Documenten/Databases/dir01-Instruments/Voice-of-dragon/03-Wav-files-predicted/04-Wav-files-calibrated-44.1kHz/modus-1-v_2filt.wav';
%       [x fs]  = Wavread(filename);
%       N = 4096; 
%       Get_LPC(x,fs,N);
% 
%       filename = '~/Documenten/Documenten-TUe/09-Training+activities/2015-Q3-Advanced-perception/College/rl001/rl001-car.wav'; 
%       [x fs]  = Wavread(filename);
%       N = 8192; 
%       Get_LPC(x,fs,N);
%
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 09/09/2014
% Last update on: 22/02/2015 % Update this date manually
% Last use on   : 22/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a section of vowel
% [x,fs]=wavread('six.wav',[24120 25930]);
% % resample to 10,000Hz (optional)

if nargin < 3
    N = 4096;
end

if fs ~= 10000
    x=resample(x,10000,fs);
    fs=10000;
end

t=(0:length(x)-1)/fs;        % times of sampling instants

% get Linear prediction filter

ncoeff=2+round(fs/1000);	% rule of thumb for formant estimation
a   = lpc(x,ncoeff);

% plot frequency response

w = Get_window('hamming',length(x));
Corr = 1;
xw = Corr*x.*w;


[hx,f]=freqz(xw,1,N,fs);
[h ,f]=freqz(1 ,a,N,fs);

x_dB = 20*log10(abs(hx)+eps);
h_dB = 20*log10(abs(h) +eps);

if nargout == 0
    
    figure;
    % plot waveform
    subplot(2,1,1);
    plot(t,x);
    legend('Waveform');
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    subplot(2,1,2);
    plot(f,h_dB); hold on
    plot(f,x_dB,'r');
    
    legend('LPC est.','x');
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
    
end

% find frequencies by root-solving

r   = roots(a);             % find roots of polynomial a
r   = r(imag(r)>0.01);   	% only look for roots >0Hz up to fs/2
ffreq=sort(atan2(imag(r),real(r))*fs/(2*pi));
                            % convert to Hz and sort

out.ncoeff = ncoeff;
out.a = a;
out.t = t;
out.f = f;

out.h_dB = h_dB;
                     
if nargout == 0
    for i=1:length(ffreq)
        fprintf('Formant %d Frequency %.1f\n',i,ffreq(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

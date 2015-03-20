function [y ydB f] = freqfft(x,K,info)
% function [y ydB f] = freqfft(x,K,info)
%
% 1. Description:
%       FFT with plot (if no output is specified)
%
%   info.typeplot - plot type (1 = semilogx, 2 = linear scale)
%       info.fs
% 
% 3. Example (note that x has to exist with a fs of 10 kHz):
%       info.fs = 10000;
%       K = 2048; % 4096-point FFT 
%       y = freqfft(x,K,info);
% 
% Programmed by Alejandro Osses, TUe
% Created on    : 08/05/2014
% Last update on: 05/06/2014
% Last use on   : 16/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 % info does not exist
    info = [];
end

if nargout == 0
    info = Ensure_field(info,'bNewFigure',1);
    info = Ensure_field(info,'typeplot',1); % 1 = semilogx
end                                         % 2 = linear scale, normalised scale
                                            % 2 = linear scale, 0 to fs/2 (if info.fs is defined)
                                            % 4 = FFT bins
% info = Ensure_field(info,'wtype',0);      % Rectangular window

if isfield(info,'fs')
    f = (1:K)/K*info.fs/2;
    if nargout == 0
        if info.typeplot == 2
            info.typeplot = 3;
        end
    end
else
    info = Ensure_field(info,'fs',44100);
    f = (1:K)/K*info.fs/2;
    warning('Default fs is considered = 44.1 kHz')
    pause(1)
end
f = f(:);

y   = fft(x,K*2);
y   = y(1:K,:);

ydB = 20*log10(abs(y));

if nargout == 0
    
    if info.bNewFigure
        figure;
    end
    switch info.typeplot
        case 1
            
            semilogx(f(:,1),ydB), grid on
            xlim([20 info.fs/2])
            xlabel(['Frequency [Hz], (fs = ' num2str(info.fs) ' [Hz])'])
        case 2
            plot((1:K)/K,ydB), grid on
            xlabel('Normalised frequency (x\pi rad/sample)')
        case 3
            plot(f,ydB), grid on
            xlabel('Frequency (Hz)')
            xlim([20 info.fs/2])
        case 4
            plot((1:K),ydB), grid on
            xlabel('FFT Bins')
    end
    ylabel('Magnitude [dB]')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
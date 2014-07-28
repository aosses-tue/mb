function y = freqfftwhpf(x,info,fc,K)
% function y = freqfftwhpf(x,info,fc,K)
%
% 1. Description:
%       FFT with plot (if no output is specified) including a High-pass 
%       filter with cut-off at fc. In case no output argument is specified,
%       an FFT is plotted
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 5/6/2014
% Last update on: 03/07/2014 % Update this date manually
% Last used on  : 03/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    K = 4096;
end

if nargin < 3
    fc = 50;
end

if nargin < 2 % info does not exist
    info = [];
end

[b,a] = hpf(fc,info.fs); % FIR filter

y = filtfilt(b,a,x);

if nargin == 0
    info = Ensure_field(info,'bNewFigure',1);
    info = Ensure_field(info,'typeplot',1);     % 1 = semilogx
                                                % 2 = linear scale, normalised scale
                                                % 2 = linear scale, 0 to fs/2 (if info.fs is defined)
% info.bNewFigure = 1;
%     freqfft(x ,K,info);
%     legend('no-HPF')
    freqfft(x1,K,info);
    legend('HPF')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
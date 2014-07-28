function y = Filters_learning(x)
% function y = Filters_learning(x)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 6/6/2014
% Last update: 6/6/2014 % Update this date manually
% Last used: 6/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    filename = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\03-Wav-files-calibrated\modus-2-v_1.wav';
end
clc

[x fs] = wavread(filename);

fc = 100;
% xlpf = 

info.fs = fs;
K = 4096;
y = freqfft(x,K,info);
fprintf('Total time RMS value: %f [dB]\n',rmsdb(x)); % equal to: 10*log10(1/ length(x)*sum(x.*x))

[yy ff] = Filterbank_analysis(x,fs,0);
% To do Filterbank: to give residuals
% To check where the cut-off freqs are...

disp(['dB values per band: ' num2str(yy')])
disp('Total energy, adding each filtered signal: ')
sum_db(yy)

% fprintf('Total time RMS value: %f [dB]\n',rmsdb(y)-10*log10(2*K)); % equal to: 10*log10(1/length(y)/K*sum(y'*y))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
function tf_speech_spectrum
% function tf_speech_spectrum
% 
% 1. Description:
%       Spectrum of normal speech according to ANSI S3.5, 1997
% 
%       Tested: yes
%       Requires that \tb_Misc\sii is added to path
% 
% Created by Tom Francart
% Edited by Alejandro Osses
% Last use: 10/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = [160 200 250 315 400 500 630 800 1000 1250 1600 2000, ...
    2500 3150 4000 5000 6300 8000];
E = SpeechSptr('normal');

figure;
plot(f,E-max(E),'b', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Relative level (dB)');
title('Spectrum of normal speech (ANSI S3.5, 1997)');
set(gca,'XScale', 'log');

% hold all
% pink = 20*log10(1./f);
% plot(f, pink-max(pink), 'r');

end
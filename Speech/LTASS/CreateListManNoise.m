function CreateListManNoise()
% function CreateListManNoise()
%
% Created by TF
% Edited by Alejandro Osses V.
% Last used on 18/06/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ListManSpectra

% 2^12 = 4096
% 2^11 = 2048
% 2^10 = 1024

%                 (order,norm-f ,amplitudes, grid to interpolate)
[H1,Fn]=freqz( fir2(2^12,F/22050,sqrt(Pxx)), 1);
[H2,Fn]=freqz( fir2(2^11,F/22050,sqrt(Pxx)), 1);
[H3,Fn]=freqz( fir2(2^10,F/22050,sqrt(Pxx)), 1);

semilogx(Fn, 10*log10(abs([H1 H2 H3])));

legend({'4096 taps' '2048 taps' '1024 taps'});

return;

%B = fir2(2^11,F/22050,sqrt(mean([Pxxnogaps2'; Pxxnogaps3'])'));
B = fir2(2^11,F/22050,sqrt(Pxx));       %% sqrt want Pxx is vermogen
w = randn(11*44100,1);
y = filter(B,1,w);
filename='/home/tom/temp/jwnoise_tom20090331.wav';
wavwrite(y,44100,16,filename);

Pxxnt = AnalyzePSDofWAV(filename,20,0);
figure;semilogx(F,10*log10(Pxxnt),F,10*log10(Pxxn),F,10*log10(Pxx));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
title('Spectra of the ListMan Lists and of the ''new'' noise');
legend({'New noise';'Speech shaped noise (jwnoise.wav)';'Speech without gaps (all lists)'});

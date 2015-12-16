clear;
clc;

[Pxxn,F] = AnalyzePSDofWAV('/data/woordenlijst/french/Bruitf.wav',20,0);
[Pxx,F] = AnalyzePSDofWAV('/data/woordenlijst/french/list/',20,0);
[Pxx2,F] = AnalyzePSDofWAV('/data/woordenlijst/french/list/',20,1e-3);

figure;
semilogx(F,20*log10(Pxx),F,20*log10(Pxxn),F,20*log10(Pxx2));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Average spectra for the French Sentences');
legend({'Speech with silent gaps (all lists)';'Speech shaped noise (bruitf.wav)';'Speech without gaps (all lists)'});

save FrenchSpectra;
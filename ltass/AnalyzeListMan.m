function AnalyzeListMan

clear;
clc;

%[Pxxn,F] = AnalyzePSDofWAV('/mnt/l/speechmaterials/dutch/jwruis.wav',20,0);
[Pxxn,F] = AnalyzePSDofWAV('/mnt/l/databases/sentences/calibration/astrid_man.wav',20,0);

[Pxx,F] = AnalyzePSDofWAV('/mnt/l/speechmaterials/dutch/LISTman/',20,0);
[Pxx2,F2] = AnalyzePSDofWAV('/mnt/l/speechmaterials/dutch/LISTman/',20,1e-3);

%[Pxx2_2048,F] = AnalyzePSDofWAV('/mnt/l/speechmaterials/dutch/LISTman/',40,1e-3);

figure;
% semilogx(F,10*log10(Pxx),F,10*log10(Pxxn),F,10*log10(Pxx2),F2,10*log10(Pxx2_2048));
semilogx(F,10*log10(Pxx),F,10*log10(Pxxn),F,10*log10(Pxx2));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Average spectra for the LIST man Sentences');
legend({'Speech with silent gaps (all lists), 4096';'Speech shaped noise (jwruis.wav)';'Speech without gaps (all lists)'; 'Speech without gaps (all lists), 2048'});

save ListManSpectra;

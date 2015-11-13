function AnalyzeNVA

pp='/mnt/l/speechmaterials/dutch/nva/all/';
[Pxx,F] = AnalyzePSDofWAV(pp,20,0,0);
[Pxx2,F2] = AnalyzePSDofWAV(pp,20,1e-3,0);

figure;
% semilogx(F,10*log10(Pxx),F,10*log10(Pxxn),F,10*log10(Pxx2),F2,10*log10(Pxx2_2048));
semilogx(F,10*log10(Pxx),F,10*log10(Pxx2));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Average spectra for the NVA CVC words');
legend({'With silent gaps';'Without silent gaps'});

save NVASpectra;


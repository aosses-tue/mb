clear;

%The 400 numbers
[Pxx,F] = AnalyzePSDofWAV('C:\temp\Numbers\',20,0,0);

%The PC format calibration noise
[Pxxn] = AnalyzePSDofWAV('C:\temp\Numbers\calibration\numbers_calibration.wav',20,0,0);

%The CD format calibration noise
[PxxnCD] = AnalyzePSDofWAV('C:\temp\Numbers\calibration\CD-track.wav',20,0,0);

figure;
hplot = semilogx(F,10*log10(Pxx),F,10*log10(Pxxn),F,10*log10(PxxnCD));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
legend({'The 400 numbers';'The calibration noise (numbers\_calibration.wav)';'The calibration noise (CD track)'});
title('Average spectrum of the numbers');
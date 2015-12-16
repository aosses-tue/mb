clear;

%The 730 sentences
[Pxx730,F] = AnalyzePSDofWAV('C:\temp\wivine\',20,0);

wavswivine;
[Pxx350,F] = AnalyzePSDofWAV(wdzfilenames,20,0);
[Pxx350_nogaps,F] = AnalyzePSDofWAV(wdzfilenames,20,1e-3);

[Pxxn] = AnalyzePSDofWAV('C:\temp\calibration\astrid_wom.wav',20,0);
[PxxnCD] = AnalyzePSDofWAV('C:\temp\calibration\CD-track-woman.wav',20,0);

figure;
hplot = semilogx(F,10*log10(Pxx730),F,10*log10(Pxx350));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
legend({'All 730 sentences';'The 350 retained sentences'});
title('Average spectrum of the sentences (gaps are retained)');

figure;
hplot = semilogx(F,10*log10(Pxx350),F,10*log10(Pxxn),F,10*log10(PxxnCD));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
legend({'The 350 sentences';'The calibration noise (astrid\_wom.wav)';'The calibration noise (CD-track)'});
title('Average spectrum of the sentences (gaps are retained)');

figure;
hplot = semilogx(F,10*log10(Pxx350_nogaps),F,10*log10(Pxx350));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
legend({'The 350 retained sentences (no gaps)';'The 350 retained sentences (gaps retained)'});
title('Average spectrum of the sentences');
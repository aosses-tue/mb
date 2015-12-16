clear;
clc;

[Pxxn,F] = AnalyzePSDofWAV('C:\temp\NVA-Children\',20,0,2);
[Pxx,F] = AnalyzePSDofWAV('C:\temp\NVA-Children\',20,0,0);


[Pxxn2,F] = AnalyzePSDofWAV('C:\temp\Gottinger I\',20,0,2);
[Pxx2,F] = AnalyzePSDofWAV('C:\temp\Gottinger I\',20,0,0);



[Pxxn3,F] = AnalyzePSDofWAV('C:\temp\Gottinger II\',20,0,2);
[Pxx3,F] = AnalyzePSDofWAV('C:\temp\Gottinger II\',20,0,0);

[Pxxnogaps,F] = AnalyzePSDofWAV('C:\temp\NVA-Children\',20,1e-3,0);
[Pxxnogaps2,F] = AnalyzePSDofWAV('C:\temp\Gottinger I\',20,1e-3,0);
[Pxxnogaps3,F] = AnalyzePSDofWAV('C:\temp\Gottinger II\',20,1e-3,0);

save ChildCDSpectra;

figure;
semilogx(F,10*log10(Pxx),F,10*log10(Pxxnogaps),F,10*log10(Pxxn));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Spectra of the LIST CD draft 2');
legend({'Speech channel (sentences) with silent gaps included';'Speech channel (sentences) without silent gaps';'Speech shaped noise sentences'});

figure;
semilogx(F,10*log10(Pxx2),F,10*log10(Pxxnogaps2),F,10*log10(Pxxn2));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Spectra of the LINT CD draft 2');
legend({'Speech channel (numbers) with silent gaps included';'Speech channel (numbers) without silent gaps';'Speech shaped noise sentences'});

figure;
semilogx(F,10*log10(Pxx3),F,10*log10(Pxxnogaps3),F,10*log10(Pxxn3));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Spectra of the LINT CD draft 2');
legend({'Speech channel (numbers) with silent gaps included';'Speech channel (numbers) without silent gaps';'Speech shaped noise sentences'});



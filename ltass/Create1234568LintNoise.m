function Create1234568LintNoise()

load LintSpectra
fs=44100;
[H1,Fn]=freqz( fir2(2^12,F/22050,sqrt(Pxx1234568)), 1, [], fs);
[H2,Fn]=freqz( fir2(2^11,F/22050,sqrt(Pxx1234568)), 1, [], fs);
[H3,Fn]=freqz( fir2(2^10,F/22050,sqrt(Pxx1234568)), 1, [], fs);

h = fdesign.arbmag('N,F,A', 2^12, F/22050, sqrt(Pxx1234568));
Hd = design(h, 'freqsamp');
[H4,Fn]=freqz( Hd.numerator, 1, [], fs);

semilogx(Fn, 20*log10(abs([H1 H2 H3 H4])), F, 10*log10(abs(Pxx1234568)));

legend({'4096 taps' '2048 taps' '1024 taps', '2048, freqsamp' 'orig' });


%B = fir2(2^11,F/22050,sqrt(mean([Pxxnogaps2'; Pxxnogaps3'])'));
B = fir2(2^11,F/22050,sqrt(Pxx1234568));       %% sqrt want Pxx is vermogen
w = randn(11*44100,1);
y = filter(B,1,w);
filename='/mnt/l/pool/tom/lint1234568tnoise_20100311.wav';
wavwrite(y,44100,16,filename);

Pxxnt = AnalyzePSDofWAV(filename,20,0);
figure;semilogx(F,10*log10(Pxxnt),F,10*log10(Pxx1234568));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
title('Spectra of the source material and of the ''new'' noise');
legend({'New noise';'Spectrum 1234568'});

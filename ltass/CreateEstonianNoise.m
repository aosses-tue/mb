function CreateEstonianNoise()

load EstonianSpectra
fs=44100;

[H1,Fn]=freqz( fir2(2^12,F/fs*2,sqrt(Pxx2)), 1);
[H2,Fn]=freqz( fir2(2^11,F/fs*2,sqrt(Pxx2)), 1);
[H3,Fn]=freqz( fir2(2^10,F/fs*2,sqrt(Pxx2)), 1);

semilogx(Fn, 10*log10(abs([H1 H2 H3])));

legend({'4096 taps' '2048 taps' '1024 taps'});

%B = fir2(2^11,F/22050,sqrt(mean([Pxxnogaps2'; Pxxnogaps3'])'));
B = fir2(2^12,F/fs*2,sqrt(Pxx2));       %% sqrt want Pxx is vermogen
w = randn(61*fs,1);
y = filter(B,1,w);

filename='/home/tom/temp/est/estnoise_20100421.wav';
wavwrite(y(fs:end),fs,16,filename);

Pxxnt = AnalyzePSDofWAV(filename,20,0);
figure;semilogx(F,10*log10(Pxx2),F,10*log10(Pxxnt));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
title('Spectra of the Est Lists and of the ''new'' noise');
legend({'Orig' 'New noise'});

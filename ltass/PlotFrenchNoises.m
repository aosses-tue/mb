function result=createFrenchNoise()

load FrenchSpectra

%B = fir2(2^11,F/22050,sqrt(mean([Pxxnogaps2'; Pxxnogaps3'])'));
% B = fir2(2^13,F/22050,sqrt(Pxx));       %% sqrt want Pxx is vermogen
% w = randn(10*44100,1);
% y = filter(B,1,w);
% filename='/data/woordenlijst/french/NewNoise.wav';
% wavwrite(y,44100,16,filename);

Pxxn=Pxxn/sqrt(mean(Pxxn.^2));
Pxx=Pxx/sqrt(mean(Pxx.^2));

%Pxxnt = AnalyzePSDofWAV(filename,20,0);
figure;semilogx(F,10*log10(Pxxn),F,10*log10(Pxx));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
title('Spectra of the French Lists and of the ''new'' noise');
legend({'Speech shaped noise (bruitf.wav)';'Speech without gaps (all lists)'});

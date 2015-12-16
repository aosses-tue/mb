function CreateAudiovizNoise()

load AudiovizSpectra

[H1,Fn]=freqz( fir2(2^12,F/fs*2,sqrt(Pxx_both)), 1);
[H2,Fn]=freqz( fir2(2^11,F/fs*2,sqrt(Pxx_both)), 1);
[H3,Fn]=freqz( fir2(2^10,F/fs*2,sqrt(Pxx_both)), 1);

semilogx(Fn, 10*log10(abs([H1 H2 H3])));

legend({'4096 taps' '2048 taps' '1024 taps'});

%B = fir2(2^11,F/22050,sqrt(mean([Pxxnogaps2'; Pxxnogaps3'])'));
B = fir2(2^11,F/fs*2,sqrt(Pxx_both));       %% sqrt want Pxx is vermogen
w = randn(61*fs,1);
y = filter(B,1,w);

% LP filter ugly distortion
[Z,P,K]=butter(15, 13500/fs*2);
sosmat=zp2sos(Z,P,K);
y=sosfilt(sosmat,y);

% video=wavread('~/dvdrip-data/audioviz/vob/001/result1+2.wav');
% rms_v=meanrms_thresh(video,44800,0.005)
%  -23.7388
y=y/meanrms(y)*10^(-23.7388/20);
filename='/home/tom/temp/audioviznoise_20100127.wav';
wavwrite(y(fs:end),fs,16,filename);

Pxxnt = AnalyzePSDofWAV(filename,20,0);
figure;semilogx(F,10*log10(Pxxnt),F,10*log10(Pxx_both));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
title('Spectra of the ListMan Lists and of the ''new'' noise');
legend({'New noise';'Speech shaped noise'});

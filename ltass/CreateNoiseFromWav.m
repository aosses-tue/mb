function CreateNoiseFromWav(sourcepath, targetfile, ntaps)

if (nargin<3)
    ntaps=2^11;
end

sourcepath=makedirend(sourcepath);
[Pxx,F] = AnalyzePSDofWAV(sourcepath,20);

figure;
hold on
semilogx(F,smooth(10*log10(Pxx)), 'b');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

B=fir2(ntaps,F/22050,sqrt(Pxx));
[H1,Fn]=freqz( B, 1, [], 44100);
semilogx(Fn, 20*log10(abs([H1])), 'r');

w = randn(11*44100,1);
y = filter(B,1,w);
wavwrite(y,44100,16,targetfile);

PxxFile = AnalyzePSDofWAV(targetfile,20,0);
semilogx(F,10*log10(PxxFile), 'g');

legend({'Overall spectrum' [num2str(ntaps) ' taps FIR filter'] 'Noise file'});
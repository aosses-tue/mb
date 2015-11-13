function AnalyzeFrenchDigitTriplet

spath='/mnt/l/speechmaterials/french/FrDigit3/broadband/';
files={};
for i=1:9
    files{end+1}=[spath num2str(i) 'a.wav'];
    files{end+1}=[spath num2str(i) 'b.wav'];
    files{end+1}=[spath num2str(i) 'c.wav'];
end

[Pxxn,F] = AnalyzePSDofWAV([spath 'frdigit3noise.wav'],20,0,2);
[Pxx,F] = AnalyzePSDofWAV(files,20,0,2);


figure;
semilogx(F,20*log10(Pxx),F,20*log10(Pxxn));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend({'FrDigit3 average', 'FrDigit3 noise'});


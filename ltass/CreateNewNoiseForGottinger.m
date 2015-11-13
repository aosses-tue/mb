clear;
clc;

load ChildCDSpectra;

%Create the New Noise signal for the Gottinger lists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = fir2(2^11,F/22050,sqrt(mean([Pxxnogaps2'; Pxxnogaps3'])'));
w = randn(5*44100,1);
y = filter(B,1,w);
wavwrite(y,44100,16,'C:\noiseI.wav');

Pxxnt = AnalyzePSDofWAV('C:\noiseI.wav',20,0);
figure;semilogx(F,10*log10(Pxxnt),F,10*log10(Pxxnogaps2),F,10*log10(Pxxnogaps3));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
title('Spectra of the Gottinger Lists (I & II) and of the ''new'' noise');
legend({'Newly designed speech shaped noise';'Gottinger I (words) without silent gaps';'Gottinger II (words) without silent gaps'});

indir = 'C:\temp\Gottinger I\';
outdir = 'C:\temp\Gottinger I new\';
d = dir([indir '*.wav']);
N = length(d);
for i=1:N
    filenames{i} = [d(i).name];
    [audio,rate] = wavread([indir filenames{i}]);
    w = randn(length(audio)+2*length(B),1);
    n = filter(B,1,w);
    clear w;
    n = n(end+1-length(audio):end);
    
    %linear gating
    n = lingate(n,0.5,0.5,rate);
    
    wavwrite([audio(:,1), n],rate,16,[outdir filenames{i}]);
end


indir = 'C:\temp\Gottinger II\';
outdir = 'C:\temp\Gottinger II new\';
d = dir([indir '*.wav']);
N = length(d);
for i=1:N
    filenames{i} = [d(i).name];
    [audio,rate] = wavread([indir filenames{i}]);
    w = randn(length(audio)+2*length(B),1);
    n = filter(B,1,w);
    clear w;
    n = n(end+1-length(audio):end);
    
    %linear gating
    n = lingate(n,0.5,0.5,rate);
    
    wavwrite([audio(:,1), n],rate,16,[outdir filenames{i}]);
end


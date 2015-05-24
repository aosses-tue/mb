function [t,s1,s2] = d_resample(inputfile,compressionFactor,opts)

if nargin < 1
    inputfile = 't53.wav';
end

if nargin < 2
    compressionFactor = 16;
end

if nargin < 3
    opts = [];
end

opts = ef(opts,'bPlot',0);
opts = ef(opts,'bSave',1);
bPlot = opts.bPlot;
bSave = opts.bSave;

[s1, sf] = wavread(inputfile);

sfNew = round(sf/compressionFactor);
t1 = (1:length(s1))/sf;
s2 = s1(1:compressionFactor:length(s1));
% sound(s2, sfNew)

if bSave
    outputfile = Delete_extension(inputfile,'wav');
    wavwrite(s2, sfNew, 16, sprintf('%s_fs_%.0f_rate_%.0f.wav',outputfile,sfNew,compressionFactor))
end
    
t2 = (1:length(s2))/sfNew;

if bPlot || bSave == 0
    figure
    subplot(2,1,1)
    plot(t1, s1)
    axis([0.1*t1(length(t1)) 0.12*t1(length(t1)) -1 1])
    subplot(2,1,2)
    plot(t2, s2)
    axis([0.1*t2(length(t2)) 0.12*t2(length(t2)) -1 1])
end

end
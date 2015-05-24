function [t,s1,s2] = d_quant(inputfile,nrOfBits,opts)

if nargin < 1
    inputfile = 't54.wav';
end

if nargin < 2
    nrOfBits = 8;
end

if nargin < 3
    opts = [];
end
opts = ef(opts,'bPlot',0);
opts = ef(opts,'bSave',1);
bPlot = opts.bPlot;
bSave = opts.bSave;

[s1, sampleFrequency] = wavread(inputfile);
t   = (1:length(s1))/sampleFrequency;
s   = round(2^(nrOfBits-1)*s1);
s2  = s/(2^(nrOfBits-1));
% sound(s2, sampleFrequency)
if bSave
    outputfile = Delete_extension(inputfile,'wav');
    wavwrite(s2, sampleFrequency, 16, sprintf('%s_Q%.0f.wav',outputfile,nrOfBits))
end

if bPlot || bSave == 0
    figure;
    subplot(2,1,1)
    plot(t, s1)
    axis([0.1*t(length(t)) 0.15*t(length(t)) -1 1])
    subplot(2,1,2)
    plot(t, s2)
    axis([0.1*t(length(t)) 0.15*t(length(t)) -1 1])
end

end
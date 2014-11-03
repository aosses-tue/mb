function [h] = Do_fluc_20141031(fi1,fi2,fi3,stPlot,fs)
% function [h] = Do_fluc_20141031(fi1,fi2,fi3,stPlot,fs)
%
% 1. Description:
%       Assesses loudness fluctuations for the audio files in fi1, fi2, fi3
% 
% 2. Stand-alone example:
%  
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file : Do_fluct_20140930
% Created on    : 31/10/2014
% Last update on: 31/10/2014 % Update this date manually
% Last use on   : 31/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    fs = 44100;
end

stPlot = Ensure_field(stPlot,'color',{'b--','r-','ko-'});
stPlot = Ensure_field(stPlot,'LineWidth',[2 1 1]);
stPlot = Ensure_field(stPlot,'Title1','L_G_{max}');

% Fluctuation strength
dBFS = 100; % Zwicker's calibration
ha = [];

if ~isnumeric(fi1)
    [insig1 fs] = Wavread(fi1);
    insig2 = Wavread(fi2);
    try
        insig3 = Wavread(fi3);
    end
else
    insig1 = fi1;
    insig2 = fi2;
    insig3 = fi3;
end

outs1 = Do_fluct(insig1,fs,dBFS);
outs2 = Do_fluct(insig2,fs,dBFS);
try
    outs3 = Do_fluct(insig3,fs,dBFS);
end

z = (1:24)-0.5;

figure;
plot( z, outs1.le_max, stPlot.color{1}, 'LineWidth',stPlot.LineWidth(1)), hold on
plot( z, outs2.le_max, stPlot.color{2}, 'LineWidth',stPlot.LineWidth(2))
try
    plot( z, outs3.le_max, stPlot.color{3}, 'LineWidth',stPlot.LineWidth(3))
end
legend(stPlot.Legend)
grid on
ylabel('Critical-band level L_G [dB]')
xlabel('Critical-band rate [Bark]')
ylim(stPlot.YLim_fig1)
title(stPlot.Title1)
% title('L_G_{max} based on N_{95}')
ha(end+1) = gca;

linkaxes(ha,'x');
xlim([0 24])

h = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [h out] = LoudnessFluctuation_offline(fi1,stPlot,fs,dBFS)
% function [h out] = LoudnessFluctuation_offline(fi1,stPlot,fs,dBFS)
%
% 1. Description:
%       Assesses loudness fluctuations for the audio files in fi1
% 
% 2. Stand-alone example:
%  
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file : Do_fluct_20140930
% Created on    : 22/01/2015
% Last update on: 16/02/2015 % Update this date manually
% Last use on   : 16/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    dBFS = 100; % Zwicker's calibration
end

if nargin < 3
    fs = 44100;
end

if nargin < 2
    stPlot = [];
end

stPlot = Ensure_field(stPlot,'color',{'b--','r-','ko-'});
stPlot = Ensure_field(stPlot,'LineWidth',[2 1 1]);
stPlot = Ensure_field(stPlot,'Title1','L_G_{max}');

ha = [];
h = [];

if ~isnumeric(fi1)
    [insig1 fs] = Wavread(fi1);
else
    insig1 = fi1;
end

outs1 = Do_fluct(insig1,fs,dBFS);

z = (1:24)-0.5;

if nargout == 0
    
    figure;
    plot( z, outs1.le_max, stPlot.color{1}, 'LineWidth',stPlot.LineWidth(1)), hold on
    plot( z, outs1.le_min, stPlot.color{2}, 'LineWidth',stPlot.LineWidth(2)), hold on
    % legend(stPlot.Legend)
    grid on
    ylabel('Critical-band level L_G [dB]')
    xlabel('Critical-band rate [Bark]')
    % ylim(stPlot.YLim_fig1)
    title(stPlot.Title1)
    % title('L_G_{max} based on N_{95}')
    ha(end+1) = gca;

    linkaxes(ha,'x');
    xlim([0 24])

    h = gcf;

end

out.t = 0;
out.z = z;
nParam = 1;
out.Data1 = outs1.le_max;
out.name{nParam} = 'Lemax';
out.param{nParam} = strrep( lower( out.name{nParam} ),' ','-');

nParam = 2;
out.Data2 = outs1.le_min;
out.name{nParam} = 'Lemin';
out.param{nParam} = strrep( lower( out.name{nParam} ),' ','-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

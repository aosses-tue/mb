function [h out outs1 outs2] = Do_fluc_20140930(fi1,fi2,stPlot,fs)
% function [h out outs1 outs2] = Do_fluc_20140930(fi1,fi2,stPlot,fs)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/09/2014
% Last update on: 24/09/2014 % Update this date manually
% Last use on   : 24/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    fs = 44100;
end

stPlot = Ensure_field(stPlot,'color',{'b--','r-'});
stPlot = Ensure_field(stPlot,'LineWidth',[2 1]);
stPlot = Ensure_field(stPlot,'Title1','L_G_{min}');
stPlot = Ensure_field(stPlot,'Title2','L_G_{max}');
stPlot = Ensure_field(stPlot,'Title3','');
stPlot = Ensure_field(stPlot,'Title4','');

% Fluctuation strength
dBFS = 100; % Zwicker's calibration
ha = [];

if ~isnumeric(fi1)
    [insig1 fs] = Wavread(fi1);
    insig2 = Wavread(fi2);
else
    insig1 = fi1;
    insig2 = fi2;
end

outs1 = Do_fluct(insig1,fs,dBFS);
outs2 = Do_fluct(insig2,fs,dBFS);

z = (1:24)-0.5;

figure;
subplot(2,1,1)
plot( z, outs1.le_min, stPlot.color{1}, 'LineWidth',stPlot.LineWidth(1)), hold on
plot( z, outs2.le_min, stPlot.color{2}, 'LineWidth',stPlot.LineWidth(2))
legend(stPlot.Legend)
grid on
ylabel('Critical-band level L_G [dB]')
ylim(stPlot.YLim_fig1)
title(stPlot.Title1)
% title('L_G_{min} based on N_{5}')
ha(end+1) = gca;

subplot(2,1,2)
plot(   z,outs1.le_min-outs2.le_min,'bo-')
ylabel('\Delta L_G [dB]')
xlabel('Critical-band rate [Bark]')
grid on
ylim(stPlot.YLim_fig2)
title(stPlot.Title3)
h = gcf;
ha(end+1) = gca;

figure;
subplot(2,1,1)
plot( z, outs1.le_max, stPlot.color{1}, 'LineWidth',stPlot.LineWidth(1)), hold on
plot( z, outs2.le_max, stPlot.color{2}, 'LineWidth',stPlot.LineWidth(2))
legend(stPlot.Legend)
grid on
ylabel('Critical-band level L_G [dB]')
ylim(stPlot.YLim_fig1)
title(stPlot.Title2)
% title('L_G_{max} based on N_{95}')
ha(end+1) = gca;

subplot(2,1,2)
plot(   z,outs1.le_max-outs2.le_max,'bo-')
ylabel('\Delta L_G [dB]')
xlabel('Critical-band rate [Bark]')
grid on
ylim(stPlot.YLim_fig2)
title(stPlot.Title4)
h(end+1) = gcf;
ha(end+1) = gca;

linkaxes(ha,'x');
xlim([0 24])

out.z       = z;
out.diff_max= outs1.le_max-outs2.le_max;
out.diff_min= outs1.le_min-outs2.le_min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

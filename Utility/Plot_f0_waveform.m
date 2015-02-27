function [h ha] = Plot_f0_waveform(t,y1,y2,options,stPlot)
% function [h ha] = Plot_f0_waveform(t,y1,y2,options,stPlot)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 26/02/2015
% Last update on: 26/02/2015 % Update this date manually
% Last use on   : 26/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h   = [];
ha  = [];

n = 2; 

stPlot = Ensure_field(stPlot,'color',{'b-','r--'});
stPlot = Ensure_field(stPlot,'LineWidth',[1 2]);

stPlot = ef(stPlot,'Title1','(a) ');
stPlot = ef(stPlot,'Title2','(b) ');
stPlot = Ensure_field(stPlot,'bYLabel',1);

figure; 
subplot(n,1,1)
plot(t,y1,stPlot.color{1}), grid on
ha = gca;

title(stPlot.Title1)
legend(options.label1);

stPlot = rmfield(stPlot,'Title1');
if stPlot.bYLabel == 1
    ylabel('Amplitude')
end

ylims = get(ha,'YLim');
set(ha(end),'YLim',1.3*ylims); % expand YLim in 20%

subplot(n,1,2) 
plot(t,y2,'r-');, grid on
ha(end+1) = gca;
title(stPlot.Title2)
legend(options.label2);

ylims = get(ha(end),'YLim');
set(ha(end),'YLim',1.3*ylims); % expand YLim in 20%

if n == 2
    xlabel('Time [s]')
end

if stPlot.bYLabel == 1
    ylabel('Amplitude')
end

h(end+1) = gcf;

linkaxes(ha,'x')
xlim(options.trange);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

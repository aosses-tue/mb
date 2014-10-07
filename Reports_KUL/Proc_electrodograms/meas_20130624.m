function meas_20130624
% function meas_20130624
%

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['Begin script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

List_of_files = {   'showsequence.m'};
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);

DirMeas     = '/home/alejandro/Documenten/Meas/';
DirThisMeas    = [DirMeas '20130624-Electrodograms/'];

addpath(DirThisMeas)
cd_old = cd;
cd(DirThisMeas)

tmax_plot = 100e-3;
t_modulator = 0:tmax_plot/10000:tmax_plot;
bPlotModulator = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1
tmin = 0.5;
ymax = 3;
putInChannel = 1;
deltaTime = 30e-3;

f = 105;
Modulator = (sin(2*pi*f*t_modulator)+1)/2;

FormatElectrodogram.Path        = DirThisMeas;
FormatElectrodogram.loop        = 'n'; % Looped recording
FormatElectrodogram.CompChannel = 16;
FormatElectrodogram.PlotGUI     = 'n';
FormatElectrodogram.Synchro     = 'n';
FormatElectrodogram.Length_analysis = 2.5;

% ACE Results:
name1 = 'SP15-104Hzm15dBFS-Tequal1-long-ACE.out';

% F0 Results:
name2 = 'SP15-104Hzm15dBFS-Tequal1-long-F0m.out';
 
name1Fig    = name2figname(name1);
name2Fig    = name2figname(name2);

[p1, p2, xx, FormatElectrodogram]    = Electrodogram_Sync(name1,name2, FormatElectrodogram);

tmax = tmin + FormatElectrodogram.time_analysis; 

h = [];

p1.current_levels = (p1.current_levels - 0)/(1) * 1;
p2.current_levels = (p2.current_levels - 0)/(1) * 1;
idx = find(p1.current_levels<100); p1.current_levels(idx) = 1;
idx = find(p2.current_levels<100); p2.current_levels(idx) = 1;

figure(1)
h(end+1) = subplot(2,1,1);
showsequence(p1), xlim([tmin tmin+tmax_plot]), ylim([0 ymax])
title(name1Fig)

h(end+1) = subplot(2,1,2);
showsequence(p2), xlim([tmin tmin+tmax_plot]), ylim([0 ymax]), hold on
title(name2Fig)

if bPlotModulator
    plot(tmin+t_modulator+deltaTime, Modulator + putInChannel, 'k')
end

linkaxes(h)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% FormatElectrodogram = [];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2

tmin = 0;
ymax = 3;
putInChannel = 1;
deltaTime = 15e-3;
f = 208;
Modulator = (sin(2*pi*f*t_modulator)+1)/2;

FormatElectrodogram.Path        = DirThisMeas;
FormatElectrodogram.loop        = 'n'; % Looped recording
FormatElectrodogram.CompChannel = 16;
FormatElectrodogram.PlotGUI     = 'n';
FormatElectrodogram.Synchro     = 'n';
FormatElectrodogram.Length_analysis = 2.5;

% ACE Results:
name1 = 'SP15-208Hzm15dBFS-Tequal1-long-ACE.out';
% 
% F0 Results:
name2 = 'SP15-208Hzm15dBFS-Tequal1-long-F0m.out';
 
name1Fig    = name2figname(name1);
name2Fig    = name2figname(name2);

[p1, p2, xx, FormatElectrodogram]    = Electrodogram_Sync(name1,name2, FormatElectrodogram);

tmax = tmin + FormatElectrodogram.time_analysis; 

h = [];

figure(2)
h(end+1) = subplot(2,1,1);
showsequence(p1), xlim([tmin tmax_plot]), ylim([0 ymax])
title(name1Fig)

h(end+1) = subplot(2,1,2);
showsequence(p2), xlim([tmin tmax_plot]), ylim([0 ymax]), hold on
title(name2Fig)

if bPlotModulator
    plot(t_modulator+deltaTime, Modulator + putInChannel, 'k')
end

linkaxes(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(cd_old);
rmpath(DirThisMeas)

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['End script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end
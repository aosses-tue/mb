function meas_20130620
% function meas_20130620
%
% Comparison between measurements done using SP12 and SP15

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['Begin script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

List_of_files = {   'showsequence.m'};
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);

DirMeas     = '/home/alejandro/Documenten/Meas/';
DirThisMeas    = [DirMeas '20130620-Physical-validation-xPC/2013-06-20-Electrodograms/'];

addpath(DirThisMeas)
cd_old = cd;
cd(DirThisMeas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1
tmin = 0;
ymax = 10;

FormatElectrodogram.Path        = DirThisMeas;
FormatElectrodogram.loop        = 'n'; % Looped recording
FormatElectrodogram.CompChannel = 16;
FormatElectrodogram.PlotGUI     = 'n';
FormatElectrodogram.Synchro     = 'n';

% ACE Results:
name1 = 'SP12-mount_R_cal.out';

% F0 Results:
name2 = 'SP15-mount_R_cal-Sens_m36.out';
 
name1Fig    = name2figname(name1);
name2Fig    = name2figname(name2);

[p1, p2, xx, FormatElectrodogram]    = Electrodogram_Sync(name1,name2, FormatElectrodogram);

tmax = tmin + FormatElectrodogram.time_analysis; 

h = [];

figure(1)
h(end+1) = subplot(2,1,1);
showsequence(p1), xlim([tmin tmax]), ylim([0 ymax])
title(name1Fig)

h(end+1) = subplot(2,1,2);
showsequence(p2), xlim([tmin tmax]), ylim([0 ymax])
title(name2Fig)

linkaxes(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FormatElectrodogram = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2

tmin = 0;
ymax = 22;

FormatElectrodogram.Path        = DirThisMeas;
FormatElectrodogram.loop        = 'n'; % Looped recording
FormatElectrodogram.CompChannel = 16;
FormatElectrodogram.PlotGUI     = 'n';
FormatElectrodogram.Synchro     = 'n';

% ACE Results:
name1 = 'SP12-mount_R_jwz-17_cut.out';

% F0 Results:
name2 = 'SP15-mount_L_jwz-17.out';
 
name1Fig    = name2figname(name1);
name2Fig    = name2figname(name2);

[p1, p2, xx, FormatElectrodogram]    = Electrodogram_Sync(name1,name2, FormatElectrodogram);

tmax = tmin + FormatElectrodogram.time_analysis; 

h = [];

figure(2)
h(end+1) = subplot(2,1,1);
showsequence(p1), xlim([tmin tmax]), ylim([0 ymax])
title(name1Fig)

h(end+1) = subplot(2,1,2);
showsequence(p2), xlim([tmin tmax]), ylim([0 ymax])
title(name2Fig)

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
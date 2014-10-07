clear, clc, close all

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['Begin script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

List_of_files = {   'showsequence.m'};
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
DirAudio = '/home/alejandro/Documenten/Meas/20130221_xPC_SP15/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference audio file:
% nameWav = 'Mix2k-ElkeZaterdag.wav';
nameWav = 'DirOut_125-250-Elke.wav';
nameAudio1 = 'VocACE_Gain40.wav';
nameAudio2 = 'VocF0m_Gain40.wav';

[x Fs]  = wavread([DirAudio, nameWav]); % 1 kHz starts at 333 ms
[x1]    = wavread([DirAudio, nameAudio1]);
[x2]    = wavread([DirAudio, nameAudio2]);
t       = 0 : 1/Fs : length(x)/Fs-1/Fs; %t = 1000*t(:); % in ms

% Reference electrodogram (Choose one and then truncate it!)

FormatElectrodogram.Path        = '/home/alejandro/Documenten/Meas/20130221_xPC_SP15/';
FormatElectrodogram.loop        = 'y'; % Looped recording
FormatElectrodogram.CompChannel = 18;
FormatElectrodogram.PlotGUI     = 'n';
FormatElectrodogram.Synchro     = 'y';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ACE Results:
% name1 = 'ACEElke-zaterdag.out';
% 
% % F0 Results:
% name2 = 'F0mElke-zaterdag.out';

% ACE Results:
name1 = 'ACE250-125-Elke.out';

% F0 Results:
name2 = 'F0m250-125-Elke.out';

name1Fig    = name2figname(name1);
name2Fig    = name2figname(name2);

[p1, p2]    = Electrodogram_Sync(name1,name2, FormatElectrodogram);

% p1.periods = 134;
% p2.periods = 139.909;

% p2.current_levels = 0.2*p2.current_levels;
tmin = 2;
tmax = tmin + 0.2;%0.2;

% figure, 
% subplot(3,1,1)
% plot(t, x), 
% xlim([0 tmax])

% subplot(1,2,1)
% showsequence(p1), xlim([tmin tmax]), ylim([0 12])
% title(name1Fig)

% subplot(1,2,2)
showsequence(p2), xlim([tmin tmax]), ylim([0 12])
title(name2Fig)

% figure,
% subplot(3,1,1)
% plot(t, x), xlim([0 tmax])
% 
% subplot(3,1,2)
% plot(t, x1), xlim([0 tmax])
% title(name1Fig)
% 
% subplot(3,1,3)
% plot(t, x2), xlim([0 tmax])
% title(name2Fig)


display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['End script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end
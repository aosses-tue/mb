clear, clc, close all

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['Begin script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

List_of_files = {   'showsequence.m'};
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
DirMeas     = '/home/alejandro/Documenten/Meas/';
DirThisMeas    = [DirMeas '20130327_xPC_Sweep/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference audio file:
% nameWav1 = 'Audio_files/Voc22-nsbF0-ACE-125+250+Elke.wav';
% nameWav2 = 'Audio_files/Voc22-nsbF0-F0m-125+250+Elke.wav';
% /home/alejandro/Documenten/Meas/20130228_xPC_SP15/

% [x1 Fs]  = wavread([DirThisMeas, nameWav1]); % 1 kHz starts at 333 ms
% [x2]     = wavread([DirThisMeas, nameWav2]);
% 
% x1 = x1(1:end-1);
% 
% t       = 0 : 1/Fs : length(x1)/Fs-1/Fs; %t = 1000*t(:); % in ms
 
% Reference electrodogram (Choose one and then truncate it!)
 
tmin = 0;
ymax = 22;

N = 1024; % N-point FFT
% f = (1:N)/N * (Fs/2);
fmin = 100;
fmax = 1e4; % 10 kHz
Hmin = -10;
Hmax = 25;

FormatElectrodogram.Path        = DirThisMeas;
FormatElectrodogram.loop        = 'n'; % Looped recording
FormatElectrodogram.CompChannel = 22;
FormatElectrodogram.PlotGUI     = 'n';
FormatElectrodogram.Synchro     = 'n';
FormatElectrodogram.time_analysis = 2; % in Seconds

tmax = tmin + FormatElectrodogram.time_analysis; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACE Results:
% name1 = 'pure125.out';
% name1 = 'sweep1000-ACE1.out';
name1 = 'sweep1seg_ACE-3.out';

% F0 Results:
% name2 = 'sweep1000-F0m3.out';
name2 = 'sweep1seg_F0mod-3.out';
 
name1Fig    = name2figname(name1);
name2Fig    = name2figname(name2);

[p1, p2]    = Electrodogram_Sync(name1,name2, FormatElectrodogram);

% p1.periods = 134;
% p2.periods = 139.909;

% % p2.current_levels = 0.2*p2.current_levels;

% 
% figure, 
% yy1 = subplot(2,1,1);
% plot(t, x1), 
% xlim([0 tmax])
% title(nameWav1)
% 
% yy2 = subplot(2,1,2);
% plot(t, x2), 
% xlim([0 tmax])
% title(nameWav2)
% 
% 
% linkaxes([yy1 yy2],'xy')

h = [];

figure
h(end+1) = subplot(2,1,1);
showsequence(p1), xlim([tmin tmax]), ylim([0 ymax])
title(name1Fig)

h(end+1) = subplot(2,1,2);
showsequence(p2), xlim([tmin tmax]), ylim([0 ymax])
title(name2Fig)

linkaxes(h)

% figure,
% h1 = freqz(x1,1,N); h2 = freqz(x2, 1, N); 
% semilogx(f, 20*log10( abs(h1) ),'b', f, 20*log10( abs(h2) ),'r'), grid on
% legend('ACE', 'F0mod')
% 
% xlim( [fmin fmax] )
% ylim( [Hmin Hmax] )

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['End script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end
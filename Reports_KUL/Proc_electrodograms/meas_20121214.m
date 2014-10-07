% Meas 2012 12 14:
%
% Dependencies:
%
%   name2figname
%   Electrodogram_Sync
%   Subplot_sequence
%   in22channel
%   AverageAmplitudePerChannel
%   FigureAverageChannel

clear, clc, close all

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('Begin script "meas_20121214.m"')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference audio file:
nameWav = 'Noises_en_delay.wav';
[x Fs]  = wavread(nameWav); % 1 kHz starts at 333 ms
t       = 0 : 1/Fs : length(x)/Fs-1/Fs; t = 1000*t(:); % in ms

% Reference electrodogram (Choose one and then truncate it!)
nameRef = 'test121214_ACE_01_trunc.out'; % 1kHz starts at 6850 ms

FormatElectrodogram.loop = 'y'; % Looped recording
FormatElectrodogram.CompChannel = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACE Results:
name1 = 'test121214_ACE_01_trunc.out';
% name2 = 'test121214_ACE_02.out';

% F0 Results:
name3 = 'test121214_F0_01.out';
% name4 = 'test121214_F0_02.out';

name1Fig    = name2figname(name1);
name3Fig    = name2figname(name3);

[p, pF0] = Electrodogram_Sync(name1,name3, FormatElectrodogram);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementing Jaime's script:

close all
Figure2 = figure; 
cant_subplot = 3; 
subplot(cant_subplot,1,1), plot(t, x), xlabel('t (ms)'), ylabel('Amplitude'), title( name2figname(nameWav) )
Subplot_sequence(p  , '', [22+1-(0:22)],0,Figure2, 2,'v',cant_subplot), title( name2figname(name1) )
Subplot_sequence(pF0, '', [22+1-(0:22)],0,Figure2, 3,'v',cant_subplot), title( name2figname(name3) )

h = ImageSetup; 
h.I_FontSize = 20; 
h.I_FontName = 'Arial'; 
h.I_Width = 8;
h.I_High= 8;
h.I_TitleInAxis = 1;
h.I_Space = [0.01,0.01];

% h.I_Ylim = [-1,1]; % Uncomment for fixing the limits in the y-axis
h.I_Xlim = [0,6000]; % 0 to 6000 ms
h.I_Grid = 'on'; 
h.I_KeepColor = 0; 
h.prepareAllFigures;

% Corriente = buffer(p.current_levels,8);


y  = in22channel(p.electrodes, p.current_levels);
y0 = in22channel(pF0.electrodes, pF0.current_levels);

Energ = AverageAmplitudePerChannel(y,6,Fs,9);
Energ0 = AverageAmplitudePerChannel(y0,6,Fs,9);
FigureAverageChannel(transpose([Energ; Energ0])/250,name1Fig,name3Fig), title('Amplitudes Normalised to 250')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subplot_sequence(y02, '', [22+1-(0:22)],0,Figure2, 2,'v',cant_subplot), title(name2Fig)

% Uncomment the following for comparing all the channels:
% % for i = 1:22
% %     [pF0, qF0, max_corr(i,1:2)] = Electrodogram_Sync_test(name4,name3, i);
% % end

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('End script "meas_20121214.m"')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
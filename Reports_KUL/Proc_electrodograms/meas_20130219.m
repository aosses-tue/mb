clear, clc, close all

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['Begin script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

List_of_files = {   'showsequence.m'};
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference audio file:
% nameWav = 'Noises_en_delay.wav';
% [x Fs]  = wavread(nameWav); % 1 kHz starts at 333 ms
% t       = 0 : 1/Fs : length(x)/Fs-1/Fs; t = 1000*t(:); % in ms

% Reference electrodogram (Choose one and then truncate it!)

FormatElectrodogram.Path        = '/home/alejandro/Documenten/Meas/20130219_Cal_xPC_SP15/';
FormatElectrodogram.loop        = 'y'; % Looped recording
FormatElectrodogram.CompChannel = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACE Results:
name1 = 'captureMinder42.out';

% F0 Results:
name2 = 'captureMinder42_F_socket_F.out';

name1Fig    = name2figname(name1);
name2Fig    = name2figname(name2);

[p1, p2]    = Electrodogram_Sync(name1,name2, FormatElectrodogram);

% p2.current_levels = 0.2*p2.current_levels;

figure, 
subplot(1,2,1)
showsequence(p1);
title(name1Fig)

subplot(1,2,2)
showsequence(p2)
title(name2Fig)

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['End script "', mfilename,'"'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end
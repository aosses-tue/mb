function [filename FS_theory_all] = ICA_FS_Create(dir_stim,bSave,dur)
% function [filename FS_theory_all] = ICA_FS_Create(dir_stim,bSave,dur)
%
% 1. Description:
%
% 2. Stand-alone example:
%       dir_stim = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-05-ICA\Stimuli\';
%       bSave = 1;
%       ICA_FS_Create(dir_stim,bSave);
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 05/05/2016
% Last update on: 05/05/2016 
% Last use on   : 14/06/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
    clc
    dir_stim = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-05-ICA\Stimuli\';
end

if nargin < 2
    bSave = 0;
end

%%%
% Common parameters:
if nargin < 3
    dur   = 2; % s
end
fs    = 44100; % Hz
N = fs*2;
FS_theory_all = [];

% 1. Reference:
f     = 1000; % Hz
fmod  = 4; % Hz
FS_theory = 1 ;
FS_theory_all = [FS_theory_all; FS_theory'];

SPL   = 60; 
mdept = 100;

filename = [];
insig = [];

for i = 1:length(fmod)
    [filename{end+1}, insig(:,end+1)] = AM_sine_savewave(f,dur,fs,fmod(i),mdept,SPL,dir_stim); 
    if bSave
        Wavwrite(insig(:,end),fs,filename{end});
    end
    % [fluct(i) fi outs] = FluctuationStrength_Garcia(insig(:,i), fs, N);
end

% 2. AM tones
fmod      = [1     2     4     8    16     32]; % Hz
FS_theory = [0.39  0.84  1.25  1.3   0.36  0.06] ;
FS_theory_all = [FS_theory_all; FS_theory'];

SPL   = 70; 

for i = 1:length(fmod)
    [filename{end+1}, insig(:,end+1)] = AM_sine_savewave(f,dur,fs,fmod(i),mdept,SPL,dir_stim); 
    if bSave
        Wavwrite(insig(:,end),fs,filename{end});
    end
    % [fluct(i) fi outs] = FluctuationStrength_Garcia(insig(:,i), fs, N);
end

% 3. AM BBN
fmod      = [1     2     4     8    16     32]; % Hz
FS_theory = [1.1205  1.5840 1.7978 1.5660 0.4793 0.1418] ;
FS_theory_all = [FS_theory_all; FS_theory'];

SPL = 60; % Hz
for i = 1:length(fmod);
    [insig(:,end+1) fname] = AM_random_noise(20,16000,SPL,dur,fs,fmod(i),mdept);
    filename{end+1} = [dir_stim fname];
    if bSave
        Wavwrite(insig(:,end),fs,filename{end});
    end
end

% 4. FM tones
f = 1500; % Hz
deltaf    = 700; % Hz
SPL       = 70; % Hz
fmod      = [1     2     4     8    16     32]; % Hz
FS_theory = [0.85  1.17  2     0.7   0.27  0.02] ;
FS_theory_all = [FS_theory_all; FS_theory'];

for i = 1:length(fmod);
    [fname, insig(:,end+1)] = FM_sine_savewave(f,dur,fs,fmod(i),deltaf,SPL,dir_stim);
    filename{end+1} = fname;
    
    if bSave
        Wavwrite(insig(:,end),fs,filename{end});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

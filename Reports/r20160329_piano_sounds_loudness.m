function r20160329_piano_sounds_loudness(dir,loud_max,notetmp)
% function r20160329_piano_sounds_loudness(dir,loud_max,notetmp)
%
% 1. Description:
%       All files inside the folder 03-Exported-as-segments will be adjusted
%       to a maximum loudness of loud_max (default = 36 sones). The new
%       audio files will be stored in a folder called '05-loudness-balanced-new'
%       Set dir to the path where the piano sounds are (see stand-alone 
%       example below). Set the target loudness loud_max to desired value.
% 
% 2. Stand-alone example:
%       dir = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\'; % close with delim
%       loud_max = 36;
%       notetmp = {'C',4};
%       r20160329_piano_sounds_loudness(dir,loud_max,notetmp);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 29/03/2016
% Last update on: 29/03/2016 
% Last use on   : 17/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bSave = 0;

if nargin == 0
    dir = [Get_TUe_data_paths('piano') '04-PAPA' delim];
end

if nargin <  2
    % 45 sone approx 95   dB
    % 36 sone approx 91.7 dB
    % 18 sone approx 81.5 dB (sine tone at 1 kHz)
    loud_max = 18; % target loudness
end
   
fstarget = 44100;
if nargin < 3
    notetmp = {'C',4}; % only C4 will be adjusted.
end

% notetmp = { 'Dsh', 1; ...
%             'F'  , 1; ... 
%             'C'  , 2; ...
%             'Ash', 2; ...
%             'F'  , 3;...
%             'C'  , 4;... 
%             'A'  , 4;...
%             'Csh', 5;... 
%             'C'  , 6;... 
%             'G'  , 6}; 

%     notetmp = { 'C'  , 2 ... 
%                 'C'  , 4, ...
%                 'A'  , 4; ... 
%                 'Csh', 5}; 

for i = 1:size(notetmp,1)

    ntmp.note   = notetmp{i,1};
    ntmp.octave = notetmp{i,2};

    noteS    = [ntmp.note num2str(ntmp.octave)];
    f0target = note2freq(ntmp);

    dirtarget{i} = sprintf('%s03-Exported-as-segments%s%s%s',dir,delim,noteS,delim); % resampled-at-44100-Hz
    gains = il_loudness_balance(dirtarget{i},fstarget,loud_max);

    if i == 1
        Mkdir([dir '05-loudness-balanced'])
    end
    dstfolder = [dir '05-loudness-balanced-new' delim noteS delim];
    Mkdir(dstfolder);

    movefile([dirtarget{i} 'loudness-balanced-new' delim],dstfolder);
    save(sprintf('%s-gains-applied4LB.mat',dstfolder),'gains');

    disp('')
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function [gains] = il_loudness_balance(directory,fs,loudness_max)
% function [gains] = il_loudness_balance(directory,fs,loudness_max)

% warning('Temporal file filter...')
file_orig = Get_filenames(directory,['*.wav']); % file_orig = Get_filenames(directory,['*.wav']);
dirout = sprintf('%sloudness-balanced%s',directory,delim);
Mkdir(dirout);

for i = 1:length(file_orig)
    filename(i) = file_orig(i);
    [insig fs] = Wavread([directory delim filename{i}]);
    
    tmp = strsplit(filename{i},'.');
    filename{i} = tmp{1}; % deletes extension...
    
    dur_ana = 350e-3; % 0.35 s
    if i == 1
        fprintf('\n\tAssumes that the piano onset occurs during the first %.1f [s] of the input signal\n\n',dur_ana);
    end
    
    gain = 0; % current gain dB
    step = 1.5; % 1.5 dB
    
    [inst_loud,short_loud] = Get_Loudness_MGB(From_dB(gain)*insig(1:round(dur_ana*fs)), fs);
    lou = max(short_loud);
    if lou > loudness_max
        go_down = 1;
    else
        go_down = 0;
    end
       
    switch go_down
        case 1
            while (lou > loudness_max & abs(lou-loudness_max) > 0.5)
                
                gain = gain - step;
                [inst_loud,short_loud] = Get_Loudness_MGB(From_dB(gain)*insig(1:round(dur_ana*fs)), fs);
                lou = max(short_loud);
                
            end
        case 0
            while ( lou < loudness_max & abs(lou-loudness_max) > 0.5 )
                
                gain = gain + step;
                [inst_loud,short_loud] = Get_Loudness_MGB(From_dB(gain)*insig(1:round(dur_ana*fs)), fs);
                lou = max(short_loud);
                
            end
        
    end
    y = From_dB(gain)*insig;
    foutname = sprintf('%s%s-%.0f-sone.wav',dirout,filename{i},loudness_max);
    Wavwrite(y,fs,foutname);
    
    gains(i) = gain;
end



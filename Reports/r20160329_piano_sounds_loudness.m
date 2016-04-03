function r20160329_piano_sounds_loudness
% function r20151127_r20160329_piano_sounds_loudness
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 29/03/2016
% Last update on: 29/03/2016 
% Last use on   : 29/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all

bSave = 0;

bDoLB           = 1; % preparing piano samples

dir = [Get_TUe_data_paths('piano') '04-PAPA' delim];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoLB
    
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
	notetmp = { 'C'  , 2 ... 
                'A'  , 4; ... 
                'Csh', 5}; 
	loud_max = [ 18, ... % 18 sones approx 81.5 dB (sine tone at 1 kHz)
                 18, ... 
                 18];
             
	fstarget = 44100;
    
    for i = 1:length(notetmp)
        
        ntmp.note   = notetmp{i,1};
        ntmp.octave = notetmp{i,2};
        
        noteS    = [ntmp.note num2str(ntmp.octave)];
        f0target = note2freq(ntmp);
                
        dirtarget{i} = sprintf('%s03-Exported-as-segments%s%s%s',dir,delim,noteS,delim); % resampled-at-44100-Hz
        gains = il_loudness_balance(dirtarget{i},fstarget,loud_max(i));
        
        if i == 1
            Mkdir([dir '05-loudness-balanced'])
        end
        dstfolder = [dir '05-loudness-balanced-new' delim noteS delim];
        Mkdir(dstfolder);
        
        movefile([dirtarget{i} 'loudness-balanced' delim],dstfolder);
        save(sprintf('%s-gains-applied4LB.mat',dstfolder),'gains');
        
    end
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function [gains] = il_loudness_balance(directory,fs,loudness_max)
% function [gains] = il_loudness_balance(directory,fs,loudness_max)

file_orig = Get_filenames(directory,['*n.wav']);
dirout = sprintf('%sloudness-balanced%s',directory,delim);
Mkdir(dirout);

for i = 1:length(file_orig)
    filename(i) = file_orig(i);
    [insig fs] = Wavread([directory delim filename{i}]);
    
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
    Wavwrite(y,fs,[dirout filename{i}]);
    
    gains(i) = gain;
end



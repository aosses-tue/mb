function r20160517_piano_sounds_loudness(dir,dir_out,loud_max)
% function r20160517_piano_sounds_loudness(dir,dir_out,loud_max)
%
% 1. Description:
%       All files inside the folder 03-Exported-as-segments will be adjusted
%       to a maximum loudness of loud_max (default = 36 sones). The new
%       audio files will be stored in a folder called '05-loudness-balanced-new'
%       Set dir to the path where the piano sounds are (see stand-alone 
%       example below). Set the target loudness loud_max to desired value.
% 
% 2. Stand-alone example:
%       dir     = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\03-Exported-as-segments\Csh5\'; % close with delim
%       dir_out = 'D:\Downloads\05-loudness-balanced\Csh5\';
%       loud_max = 18;
%       r20160517_piano_sounds_loudness(dir,dir_out,loud_max);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 29/03/2016
% Last update on: 29/03/2016 
% Last use on   : 17/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    dir = uigetdir(pwd);
    switch dir
            case 0
                error('Please select a valid directory. Re-run this code...')
        otherwise % then a directory has been found
            dir = [dir delim];
    end
end

if nargin < 2
    dir_out = uigetdir(pwd);
    switch dir_out
            case 0
                error('Please select a valid directory where to write the loudness balanced stimuli. Re-run this code...')
        otherwise % then a directory has been found
            dir_out = [dir_out delim];
    end
end

if nargin <  3
    % 45 sone approx 95   dB
    % 36 sone approx 91.7 dB
    % 18 sone approx 81.5 dB (sine tone at 1 kHz)
    loud_max = input('Please enter the target loudness value (default of 18 sones): '); % 18; % target loudness
end
   
gains = il_loudness_balance(dir,dir_out,loud_max);
save(sprintf('%s-gains-applied4LB.mat',dir_out),'gains');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function [gains] = il_loudness_balance(dir,dir_out,loudness_max)
% function [gains] = il_loudness_balance(directory,fs,loudness_max)

file_orig = Get_filenames(dir,['*.wav']); % file_orig = Get_filenames(directory,['*.wav']);

for i = 1:length(file_orig)
    filename(i) = file_orig(i);
    [insig fs] = Wavread([dir delim filename{i}]);
    
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
    foutname = sprintf('%s%s-%.0f-sone.wav',dir_out,filename{i},loudness_max);
    Wavwrite(y,fs,foutname);
    
    gains(i) = gain;
end



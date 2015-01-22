function VoD_convert_txt2wav
% function VoD_convert_txt2wav
%
% 1. Summary:
%       Converts txt-files containing 3-channel 'voice of the dragon' 
%       measurements into wav-files.
%       Note that *-3 files always saturate, however, these signals were 
%       recorded just to determine period of rotation
% 
% 2. Stand-alone example:
%       VoD_convert_txt2wav;
%
% 3. Additional info:
%   Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/05/2014
% Last update on: 28/07/2014 % Update this date manually
% Last use on   : 28/07/2014 % Update this date manually
% 
% Original file name: Convert_txt2wav, changed on: 27/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;

if bDiary
    Diary(mfilename)
end

main_dir    = Get_TUe_paths('db_voice_of_dragon');
output_dir  = [main_dir '02-Wav-files' delim 'new' delim];
fs          = 10000;
Att         = 1; % From_dB(-20);

dirs = {['1 referentie'     delim], ...
        ['2 omgedraaid'     delim], ...
        ['3 schuim hand'    delim], ...
        ['4 besneden'       delim], ...
        ['5 hornetje er af' delim], ...
        ['6 ledenmaat eraf' delim]};

for i = 1:length(dirs)
    
    Mkdir([output_dir dirs{i}]); % Creates output directory in case it does not exist
    
    files = dir([main_dir dirs{i} '*.txt']);
    
    for j = 1:length(files)
        
        y = import_physical_measure([main_dir dirs{i} files(j).name], 1);
        
        Wavwrite(Att*y(:,1),fs, [output_dir dirs{i} Delete_extension(files(j).name,'txt') '-1.wav']);
        Wavwrite(Att*y(:,2),fs, [output_dir dirs{i} Delete_extension(files(j).name,'txt') '-2.wav']);
        Wavwrite(Att*y(:,3),fs, [output_dir dirs{i} Delete_extension(files(j).name,'txt') '-3.wav']);
        
    end
    
end

if bDiary
    diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
function VoD_one_predicted_file_from_txt(filename)
% function VoD_one_predicted_file_from_txt(filename)
%
% 1. Description:
%       Generate predicted Voice of the Dragon audio files. The data were 
%       previously generated in txt files from Mathematica
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3.1 Stand-alone example to generate all wav-files:
%       VoD_predicted_files_from_txt;
%
% 3.2 Stand-alone example to generate Mode 2, near-field:
%       do_modes    = 2;
%       do_fields   = 2;
%       VoD_predicted_files_from_txt(do_modes, do_fields);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 27/05/2014
% Last update on: 18/06/2014 % Update this date manually
% Last use on   : 18/06/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;

if bDiary
    Diary(mfilename)
end

info.bSave  = input([mfilename '.m: type 1 if you want to generate VoD data: ']);
paths.VoD   = Get_TUe_paths('db_voice_of_dragon');
bDoCheckCompatibility = 0;

dir_initial = [paths.VoD '03-Wav-files-predicted' delim '01-Model' delim 'Data' delim];

if nargin < 1
    [f1 f2] = uigetfile(dir_initial);
    filename = [f2 f1];
    dir_predicted = f2;
end

fs      = 10000; % known from Mathematica model
dBSPL   = 65;

if bDoCheckCompatibility
    
    disp('Are you sure you want to check compatibility?');
    disp('This is going to take about 30-35 minutes per file');
    disp('10 files = approx. 5 to 6 Hours')
    disp('Press any button to continue...')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause;
    
    Check_txt_compatibility(filename);
	% compatible txt files are stored with suffix '-c'
    
end
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% filename = [dir_predicted 'mode-' num2str(k-1) '-v_' num2str(fieldtype) '-c.txt'];
outfilename = Delete_extension(filename,'txt');
try
    x = import_physical_measure(filename, 1,inf,1);
catch
    error(['Check if all the (corrected files) wav-files are stored at ' dir_predicted '.\nMost probably you have to set ''bDoCheckCompatibility to 1'''])
end

tx      = (1:length(x))/fs;
y       = setdbspl(x,dBSPL);

if info.bSave
    if fs ~= 44100
        y = resample(y,44100,fs);
        fs = 44100;
    end
    Wavwrite(y,fs,outfilename)
else
    disp([outfilename ' loaded but not stored'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDiary
    diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
function VoD_predicted_files_from_txt(do_modes, do_fields)
% function VoD_predicted_files_from_txt(do_modes, do_fields)
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

bDiary = 1;

if bDiary
    Diary(mfilename)
end

info.bSave  = input([mfilename '.m: type 1 if you want to generate VoD data: ']);
paths.VoD   = Get_TUe_paths('db_voice_of_dragon');
bDoCheckCompatibility = 1;

if nargin < 1
    do_modes    = 2:5;
end

if nargin < 2
    do_fields   = 1:2;
end

dir_initial = [paths.VoD '03-Wav-files-predicted' delim '01-Model' delim 'Data' delim];

dir_predicted = uigetdir(dir_initial,'Select the directory where the Mathematica''s txt files are stored: ');
dir_predicted = [dir_predicted delim];

fs      = 10000; % known from Mathematica model
dBSPL   = 65;

if bDoCheckCompatibility
    
    disp('Are you sure you want to check compatibility?');
    disp('This is going to take about 30-35 minutes per file');
    disp('10 files = approx. 5 to 6 Hours')
    disp('Press any button to continue...')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause;
    
    for fieldtype = do_fields
        for k = do_modes % Mode from 2 to 5

            filename = [dir_predicted 'mode-' num2str(k-1) '-v_' num2str(fieldtype) '.txt'];
            Check_txt_compatibility(filename);
            % compatible txt files are stored with suffix '-c'
        end
    end
end
              
for fieldtype = do_fields % 1 = far-field, 2 = near-field 
    for k = do_modes % Mode from 2 to 5

        filename = [dir_predicted 'mode-' num2str(k-1) '-v_' num2str(fieldtype) '-c.txt'];
        outfilename = [dir_predicted 'modus-' num2str(k-1) '-v_' num2str(fieldtype)];
        try
            x = import_physical_measure(filename, 1,inf,1);
        catch
            error(['Check if all the (corrected files) wav-files are stored at ' dir_predicted '.\nMost probably you have to set ''bDoCheckCompatibility to 1'''])
        end

        tx      = (1:length(x))/fs;
        y       = setdbspl(x,dBSPL);

        if info.bSave
            Wavwrite(y,fs,outfilename)
        else
            disp([outfilename ' loaded but not stored'])
        end
    end
end

if info.bSave == 1
    bHPF = 1;
    options.bSave = 1;
    VoD_run(bHPF,options,[],do_modes); %,info,take,modes2check) % To make sure that calibrated audio files are generated
end

if bDiary
    diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
function VoD_generate_formants(bPraat)
% function VoD_generate_formants(bPraat)
%
% 1. Description:
%       If bPraat is 1 then new formants from Praat are going to be assessed
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       VoD_generate_formants(bPraat);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 12/08/2014
% Last update on: 12/08/2014 % Update this date manually
% Last use on   : 12/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bPraat = 1;
end

% paths = Get_TUe_paths;
subpaths = Get_TUe_subpaths('db_voice_of_dragon');

% src_measured = subpaths.dir_calibrated_m;
% src_modelled = subpaths.dir_calibrated_p;
src_measured = subpaths.dir_calibrated_ms;
src_modelled = subpaths.dir_calibrated_ps;

if bPraat
    % Get formants for measured
    Get_formants_from_dir(src_measured);

    % Get formants for modelled
    Get_formants_from_dir(src_modelled);
end

% Move files:
status1 = Mkdir(subpaths.dir_fn_m);
status2 = Mkdir(subpaths.dir_fn_p);

if ~status1 | ~status2
    disp('You might be overwriting files, do you want to proceed, press any key (ctrl+c to cancel)...')
    pause();
end

% file1log    = Get_filenames(subpaths.dir_calibrated_m,'log*');
% files1      = Get_filenames(subpaths.dir_calibrated_m,'*filt.txt');
% file2log    = Get_filenames(subpaths.dir_calibrated_p,'log*');
% files2      = Get_filenames(subpaths.dir_calibrated_p,'*filt.txt');

file1log    = Get_filenames(src_measured,'log*');
files1      = Get_filenames(src_measured,'*filt.txt');
file2log    = Get_filenames(src_modelled,'log*');
files2      = Get_filenames(src_modelled,'*filt.txt');

for i = 1:length(files1)
    movefile([src_measured files1{i}], [subpaths.dir_fn_m files1{i}]);
end

for i = 1:length(files2)
    movefile([src_modelled files2{i}], [subpaths.dir_fn_p files2{i}]);
end
movefile([src_measured file1log{1}], [subpaths.dir_fn_m file1log{1}]);
movefile([src_modelled file2log{1}], [subpaths.dir_fn_p file2log{1}]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

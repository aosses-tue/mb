function info = Get_F0_AC_praat_from_dir(directory, bDiary)
% function info = Get_F0_AC_praat_from_dir(directory, bDiary)
%
%  1. Description:
%       Creates F0 contours using Praat (in case they do not exist yet) for 
%       all the audio files (*.wav) under the info.wav directory.
%
%   2. Additional info:
%       
%   3. Examples:
%   3.1 Old example (valid for versions before 31/07/2014):
%       info = analysis_F0_wav;
% 
%   3.2 Example. Run the following line, then select the directory:
%       info = Get_F0_AC_praat_from_dir;
% 
%   3.3 Example
%       dirs = Get_TUe_subpaths('db_voice_of_dragon');
%       Get_F0_AC_praat_from_dir(dirs.dir_calibrated_m); % F0 Extraction of measured signals
%       Get_F0_AC_praat_from_dir(dirs.dir_calibrated_p); % F0 Extraction of predicted signals
%    
% Programmed by Alejandro Osses, ExpORL, Belgium, 2014
% Original name : analysis_F0_wav (changed on 31/07/2014)
% Created on    : March 2014
% Last update on: 31/07/2014
% Last use on   : 31/07/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2 
    bDiary = 1;
end

info = [];
if nargin == 0
    directory = [uigetdir delim];
end

if nargin < 2
    p = [];
end

p = Ensure_field(p, 'F0max', 1400); % default KU Leuven = 300 Hz
p = Ensure_field(p, 'F0min',   75);

close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info    = Ensure_field(info, 'speaker', ''); % VlMatrix 
info    = Ensure_field(info, 'tstep'  , 10e-3);
info    = Ensure_field(info, 'analysed_wav_files',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info            = Ensure_field(info, 'root_dir'  , directory);

extra.bExtension = 0; % To delete extension
file_orig       = Get_filenames(info.root_dir,['*.wav'],extra);

if length(file_orig)~=0
    
    try
        tmp = Get_date;
        info.outputs  = [Get_TUe_paths('outputs') tmp.date2print delim];
    catch
        info.outputs  = [info.root_dir 'F0_from_praat'    delim]; 
    end
    Mkdir(info.outputs);
    Diary(mfilename,bDiary,info.outputs);
    
    fprintf('\n\n');
    fprintf('Running script %s.m\n', mfilename);
    fprintf('\tF0 Extraction of files inside dir %s (input 1)\n', info.root_dir);
    fprintf('\tResults are going to be stored at %s (see log-file)\n\n', info.outputs);
    fprintf('List of files to be processed:\n')
    fprintf('\t#\t file name\n')
    for k = 1:length(file_orig)
        fprintf('\t%.0f\t%s.wav\n',k,file_orig{k});
    end
    
end

close all, clc
cumF0   = 0;
cumtF0  = 0;
cumt_uv = 0;

for k = 1:length(file_orig)
    
    disp(['Processing: ' file_orig{k}])
    
    [file bCompatible] = Has_filename_without_spaces(file_orig{k});

    if bCompatible == 1
        filename    = [info.root_dir delim file];
    else
        disp([mfilename '.m: a non-compatible file name was found...copying and re-naming files...'])
        dir_wav_backup  = [info.root_dir delim 'wav' delim];
        filename_i      = [info.root_dir delim file_orig{k}];
        filename        = [dir_wav_backup      file];
        Mkdir( dir_wav_backup )
        copyfile([filename_i '.wav'], [filename '.wav'])
    end

    info.analysed_wav_files = [info.analysed_wav_files; {file}];

    [y fs]      = wavread([filename '.wav']); % Clean audio file, reference

    t           = (1:length(y))/fs; % Only used to plot

    %%% Loading F0 contour from Praat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    outputfile  = [info.outputs file '.txt'];

    try
        [tF0 F0ref] = Get_F0_praat_from_txt(outputfile,p);

    catch % Then txt files from Praat do not exist yet

        inputfile = [filename '.wav'];
        type = 3; 
        [tF0 F0ref]    = Get_F0_AC_praat(inputfile,outputfile,type);
        if type == 0 | type == 1 | type == 2
            disp('Not using TU/e approach')
        end

    end

    %%% end of 'Loading F0 contour from Praat' %%%%%%%%%%%%%%%%%%%%

    idx         = find(F0ref < p.F0min);
    F0ref(idx)  = NaN;

    tref        = 0 + info.tstep: info.tstep : max(t); %max(tF0); % Same time-span than *.fx, but with fixed time steps
    F0ref       = interp1(tF0,F0ref, tref); % F0 is resampled
    idx         = find(F0ref > p.F0max);
    F0ref(idx)  = NaN;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Database info:
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cum_idx = find(F0ref > 0);

    cumtF0  = cumtF0 + sum(F0ref > 0);
    cumF0   = cumF0 + sum(F0ref(cum_idx));
    cumt_uv = cumt_uv + sum(isnan(F0ref));
    Avg_F0_tot  = cumF0/cumtF0;
    Avg_F0 = mean(F0ref(cum_idx));
    local_F0min = min(F0ref);
    local_F0max = max(F0ref);

end

% Displays some information:

disp(['(F0min, F0max)[Hz]: (' num2str(p.F0min) ', ' num2str(p.F0max) ')']);
disp([directory ' database, ' info.speaker, ': '])
disp(['Total time     [s]: ' num2str(info.tstep*(cumtF0+cumt_uv))])
disp(['Total voiced   [s]: ' num2str(info.tstep*(cumtF0))])
disp(['Total unvoiced [s]: ' num2str(info.tstep*(cumt_uv))])
disp(['Avg F0        [Hz]: ' num2str(Avg_F0_tot)])
disp(['Number of analysed audio files: ' num2str(k)])

if bDiary
    diary off
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

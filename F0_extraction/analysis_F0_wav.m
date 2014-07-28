function info = analysis_F0_wav(directory)
% function info = analysis_F0_wav(directory)
%
% Creates F0 contours using Praat (in case they do not exist yet) for all 
% the audio files (*.wav) under the info.wav directory
%
% % Example:
%       info = analysis_F0_wav;
% 
% Programmed by Alejandro Osses, ExpORL, 2014
% Created on     : March 2014
% Last updated on: 19/05/2014
% Last used on   : 19/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = [];
if nargin == 0
    directory = [uigetdir delim];
end

if nargin < 2
    p = [];
end

p = Ensure_field(p, 'F0max', 300);
p = Ensure_field(p, 'F0min',  75);

close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% info    = Ensure_field(info, 'speaker', 'LIST-f'); % Female 
info    = Ensure_field(info, 'speaker', ''); % VlMatrix 
info    = Ensure_field(info, 'tstep'  , 10e-3);
info    = Ensure_field(info, 'analysed_wav_files',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info            = Ensure_field(info, 'root_dir'  , directory);
% info            = Ensure_field(info, 'results_dir', [info.root_dir 'results' delim]);

info.praat      = [info.root_dir 'praat'    delim]; 

extra.bExtension = 0; % To delete extension
file_orig       = Get_filenames(info.root_dir,['*.wav'],extra);

if length(file_orig)~=0
    Mkdir(info.praat);
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

    outputfile  = [info.praat file '.txt'];

    try
        [tF0 F0ref] = Get_F0_praat_from_txt(outputfile);

    catch % Then txt files from Praat do not exist yet

        inputfile = [filename '.wav'];
        type = 0; 
        [tF0 F0ref]    = Get_F0_AC_praat(inputfile,outputfile,type);
        if type == 1 | type == 2
            disp('Using de Cheveigne''s aprroach')
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
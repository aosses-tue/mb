function misc = Get_TUe_subpaths(type) 
% function misc = Get_TUe_subpaths(type)
%   
%   Get local subpaths.
% 
%   type:
%       'db_voice_of_dragon'
%       'tb_Loudness_v12'
%       'tb_NMT' (Nucleus MATLAB Toolbox)
%   
%   Where:
%       db = database
%       ex = experiment
%       lx = LaTeX
%       tb = toolbox
%
% % Example 1: Subpaths under NMT toolbox:
%       misc_sub = Get_TUe_subpaths('tb_NMT');
%
% Programmed by Alejandro Osses, TUe, 2014
% Created on    : 24/06/2014
% Last update on: 31/07/2014
% Last use on   : 31/07/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = Get_TUe_paths(type);

if strcmp(type,'db_voice_of_dragon')
    
    misc.dir_meas_wav      = [path '02-Wav-files' delim];
    misc.dir_calibrated_m  = [misc.dir_meas_wav '03-Wav-files-calibrated' delim]; 
    misc.dir_calibrated_ms = [misc.dir_meas_wav '04-Wav-files-calibrated-synchro' delim];
    misc.dir_predicted_txt = [path '03-Wav-files-predicted' delim '01-model' delim 'Data' delim];
    misc.dir_calibrated_p  = [path '03-Wav-files-predicted' delim '03-Wav-files-calibrated' delim]; 
    misc.dir_calibrated_ps = [path '03-Wav-files-predicted' delim '05-Wav-files-calibrated-synchro' delim]; 
    misc.dir_f0_m          = [path '04-fn-extraction' delim 'measured_f0' delim];
    misc.dir_f0_p          = [path '04-fn-extraction' delim 'modelled_f0' delim];
    misc.dir_fn_m          = [path '04-fn-extraction' delim 'measured_fn' delim];
    misc.dir_fn_p          = [path '04-fn-extraction' delim 'modelled_fn' delim];
    
    misc.dir_measurements  ={[misc.dir_meas_wav '1 referentie'     delim], ...
                             [misc.dir_meas_wav '2 omgedraaid'     delim], ...
                             [misc.dir_meas_wav '3 schuim hand'    delim], ...
                             [misc.dir_meas_wav '4 besneden'       delim], ...
                             [misc.dir_meas_wav '5 hornetje er af' delim], ...
                             [misc.dir_meas_wav '6 ledenmaat eraf' delim]};
                             
    misc.dir_meas_def      = misc.dir_measurements{1}; % Anechoic condition

elseif strcmp(type,'tb_Loudness_v12')
    
    misc.WAV               = [path 'WAV' delim];
    
elseif strcmp(type,'tb_NMT')
    
    misc.Processing        = [path 'Matlab' delim 'Processing' delim]; % Ensure_field is here
    misc.Sequence          = [path 'Matlab' delim 'Sequence'   delim]; % Collate_into_sequence.m
    misc.Utility           = [path 'Matlab' delim 'Utility'    delim]; % From_dB; To_dB; Matlab_version.m
    
elseif strcmp(type,'tb_Plot4papers_JU')
    misc.Utility           = [path 'Utility'        delim];
else
        misc = 4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
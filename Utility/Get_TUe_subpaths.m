function misc = Get_TUe_subpaths(type) 
% function misc = Get_TUe_subpaths(type)
%   
% 1. Description:
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
% Last update on: 17/09/2014
% Last use on   : 20/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = Get_TUe_paths(type);

if strcmp(type,'db_speechmaterials')
    
    misc.LISTf  = [path 'dutch' delim 'list' delim]; 
    misc.allfiles_PB = [path 'english' delim 'paul-bagshaw' delim 'wav'      delim];
    misc.fda_eval_PB = [path 'english' delim 'paul-bagshaw' delim 'fda_eval' delim];
    misc.allfiles_LISTf = [path 'dutch' delim 'list' delim 'alle-zinnen' delim];
    misc.fda_eval_LISTf = [path 'dutch' delim 'list' delim 'fda_eval' delim];
    misc.fda_eval_LISTf_wav_info = [misc.fda_eval_LISTf 'wav_info'];

elseif strcmp(type,'db_voice_of_dragon')
    
    misc.dir_meas_wav      = [path '02-Wav-files' delim];
    misc.dir_calibrated_m  = [misc.dir_meas_wav '03-Wav-files-1-referentie'         delim]; % before 20/01/2015 = '03-Wav-files-calibrated'
    misc.dir_calibrated_m6 = [misc.dir_meas_wav '03-Wav-files-6-lederenmaat-eraf'   delim]; % created on 20/01/2015
    
    misc.dir_calibrated_ms = [misc.dir_meas_wav '04-Wav-files-calibrated-synchro'   delim];
    misc.dir_predicted_txt = [path '03-Wav-files-predicted' delim '01-model' delim 'Data' delim];
    misc.dir_calibrated_p  = [path '03-Wav-files-predicted' delim '03-Wav-files-1-referentie' delim]; % before 20/01/2015 = '03-Wav-files-calibrated'
    misc.dir_calibrated_p6 = [path '03-Wav-files-predicted' delim '03-Wav-files-6-lederenmaat-eraf' delim]; % created on 20/01/2015
    
    misc.dir_calibrated_ps = [path '03-Wav-files-predicted' delim '05-Wav-files-calibrated-synchro' delim]; 
    misc.dir_wav_all       = [path '04-Wav-files-all' delim];
    misc.dir_f0_m          = [path '05-fn-extraction' delim 'measured_f0' delim];
    misc.dir_f0_p          = [path '05-fn-extraction' delim 'modelled_f0' delim];
    misc.dir_fn_m          = [path '05-fn-extraction' delim 'measured_fn' delim];
    misc.dir_fn_p          = [path '05-fn-extraction' delim 'modelled_fn' delim];
    
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
    
elseif strcmp(type,'Reports_KUL')
    
    misc.Experiments = [path 'Experiments' delim];
    misc.Proc_electrodograms = [path 'Proc_electrodograms' delim];
        
else
        misc = 4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
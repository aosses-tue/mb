function misc = Get_TUe_paths(type) 
% function misc = Get_TUe_paths(type)
%   
%   Get local paths as structured at the TU/e Dell laptop (since 01/05/2014).
%   If 'type' is specified, then 'misc' will be a character containing the
%   that specific path. See example 2 (below).
% 
%   types available:
%       'db_fastl2007'
%       'db_voice_of_dragon'
%       'db_speechmaterials' - location X-drive
%       'db_speechmaterials_local' - local back-up of speech materials
%       'db_voice_of_dragon'
%       'Music_acoustics'
%       'praat'     - location of Praat (terminal: whereis praat)
%       'tb_AM'     - Auditory Modelling Toolbox
%       'tb_NMT'    - Nucleus MATLAB Toolbox    
%   
%   Where:
%       db = database
%       ex = experiment
%       lx = LaTeX
%       tb = toolbox
% 
% % Example 1: to get all important directories
%       misc = Get_TUe_paths;
%
% % Example 2: to get directory where database of the voice of the dragon is:
%       misc = Get_TUe_paths('db_voice_of_dragon');
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 15/03/2014
% Last update on: 07/02/2016
% Last use on   : 10/02/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alejandro's paths

info.bOnline = 0;
    
if nargout == 0 & nargin == 0
    disp([mfilename '.m= Change bOnline to 0 in case you don''t have success if a directory is not found'])
    pause(2)
end

if isunix
    
    misc.username           = 'alejandro';
    misc.MATLAB_KUL         = ['~/Documenten/MATLAB/MATLAB_svn/']; % KU Leuven MATLAB
    
    if info.bOnline
        misc.MATLAB         = ['~/Documenten/Dropbox/TUe/MATLAB/']; % TU/e MATLAB dropbox
    else
        % misc.MATLAB         = ['~/Documenten/MATLAB/MATLAB_TUe/']; % commented on 03/08/2014
        misc.MATLAB         = ['~/Documenten/MATLAB/MATLAB_git/'];
    end
    misc.Databases          = ['~/Documenten/Databases' delim];

    misc.outputs            = ['~/Documenten/MATLAB/outputs' delim];
    misc.praat              =  '/usr/bin/praat'; % whereis praat
    misc.praat_scripts      = ['~/Documenten/Praat_svn/'];
    misc.SVN_KUL            = misc.MATLAB_KUL;
    % % Other KU Leuven paths (add them if necessary):
    % misc.svn.Meas           = ['/home/' misc.username '/Documenten/Meas/Meas/'];
    % misc.result_folder      = [misc.svn.Meas 'Experiments/Results_XML/'];
    % misc.experiment_report  = [misc.svn.Meas 'Experiments/Experiment_reports/'];
    % misc.Stats              = [misc.MATLAB_KUL 'Statistics'   delim];
else
    
    % % WinXP, VirtualBox Alejandro on 28.08.2013:
    % misc.username   =  'Administrator';
    % misc.MATLAB     = ['C:\Documents and Settings\' misc.username '\Desktop\alejandro\MATLAB_svn' delim];
    
    % % Win7, TU/e:
    misc.username           =  'aosses';
    if info.bOnline
        misc.MATLAB             = ['G:\MATLAB' delim];
    else
        % misc.MATLAB         = ['D:\MATLAB-off-line' delim];
        misc.MATLAB         = ['D:\MATLAB_git' delim];
    end
    misc.MATLAB_KUL         = ['D:\SVN-KU-Leuven\alejandro\MATLAB_svn' delim];
    misc.Databases          = ['D:\Databases' delim];
    misc.SVN_KUL            = ['D:\SVN-KU-Leuven\alejandro' delim];
    
    misc.CIs                =  [misc.MATLAB 'CIs' delim];
    
    misc.ex_APEX_results    = [misc.SVN_KUL 'Meas' delim 'Experiments' delim 'Results_XML' delim]; % work at KUL
    misc.outputs            = ['D:\Output' delim];
    misc.praat              = ['C:\praat5376_win32\praatcon.exe'];
    
end

misc.DSP                = [misc.MATLAB   'DSP'             delim]; 
misc.F0_extraction      = [misc.MATLAB   'F0_extraction' delim];
misc.Filterbank         = [misc.MATLAB   'Filterbank'      delim];
misc.language           = 'NL'; % other possibilities: 'EN'
misc.Localisation       = [misc.MATLAB   'Localisation'    delim];
misc.praat_scripts      = [misc.MATLAB   'Praat'           delim];
misc.PA_stim_generation = [misc.MATLAB   'PA_stim_generation' delim];
misc.Psychoacoustics    = [misc.MATLAB   'Psychoacoustics' delim];
misc.FluctuationStrength= [misc.Psychoacoustics 'FluctuationStrength_Garcia' delim];
misc.FluctuationStrength_Sontacchi= [misc.Psychoacoustics 'FluctuationStrength_Sontacchi' delim];
misc.Reports            = [misc.MATLAB   'Reports'         delim];
misc.Reports_KUL        = [misc.MATLAB   'Reports_KUL'     delim];
misc.Roughness_Daniel   = [misc.Psychoacoustics 'Roughness_Daniel'   delim];
misc.Roughness_Duisters = [misc.Psychoacoustics 'Roughness_Duisters' delim];
misc.Speech             = [misc.MATLAB   'Speech'          delim];
misc.ltass              = [misc.Speech   'LTASS'           delim];
misc.ICRA               = [misc.Speech   'ICRA'            delim];
misc.ICRA_Tobias        = [misc.Speech   'ICRA_Tobias'     delim];
misc.Misc               = [misc.MATLAB   'tb_Misc'         delim];
misc.Music_acoustics    = [misc.MATLAB   'Music_acoustics' delim];
misc.pl_Classes         = [misc.MATLAB_KUL 'Classes'       delim]; % Matthias' classes to plot
misc.sii                = [misc.Misc     'sii'             delim];
misc.Simulink           = [misc.MATLAB_KUL 'Simulink_XPC'  delim];
misc.Stats              = [misc.MATLAB   'Stats'           delim];
misc.Text               = [misc.MATLAB  'Text'             delim];
misc.tb_AM              = [misc.MATLAB   'tb_AM'           delim];
misc.tb_AM_AddOns       = [misc.MATLAB   'tb_AM_AddOns'    delim];
misc.tb_AM_my_examples  = [misc.tb_AM    'my_examples'     delim]; % Added on 07/10/2014
misc.tb_Misc            = [misc.MATLAB   'tb_Misc'         delim];
misc.tb_NMT             = [misc.MATLAB   'tb_NMT_4.31'     delim];
misc.tb_NMTAddOns       = [misc.MATLAB   'tb_NMT_AddOns'   delim]; 
misc.tb_Plot4papers_JU  = [misc.MATLAB   'tb_Plot4papers_JU' delim]; % Developed by Jaime Undurraga
misc.tb_Psysound32      = [misc.MATLAB   'tb_Psysound32'   delim];
misc.tb_SP15            = [misc.MATLAB   'tb_SP15'         delim];
misc.tb_SP15AddOns      = [misc.MATLAB   'tb_SP15_AddOns'  delim];
% misc.FluctuationStrength_TestTones = [misc.outputs 'FS_test_tones' delim]; % temporal folder

misc.tb_APEX            = [misc.MATLAB   'tb_APEX'         delim];
misc.tb_APEX_templates  = [misc.tb_APEX  'a3templates'     delim];
misc.tb_APEX_AddOns     = [misc.MATLAB   'tb_APEX_AddOns'  delim];
misc.tb_APEX_tools      = [misc.tb_APEX  'tools'           delim];
misc.tb_Loudness_v12    = [misc.MATLAB   'tb_Loudness_v12' delim];

misc.Utility_KUL        = [misc.MATLAB   'Utility_KUL'     delim];
     
if (nargin==1)
    if (~isfield(misc,type))
        error('Invalid type');
    end
    misc=misc.(type);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

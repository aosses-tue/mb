function misc = Get_TUe_data_paths(type)
% function misc = Get_TUe_data_paths(type)
%
% 1. Description:
%
% 2. Stand-alone example:
%       dir_audacity = [Get_TUe_data_paths('db_audacity') delim];
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/08/2014
% Last update on: 28/10/2014 % Update this date manually
% Last use on   : 28/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alejandro's paths

info.bOnline = 0;
    
if nargout == 0 & nargin == 0
    disp([mfilename '.m= Change bOnline to 0 in case you don''t have success if a directory is not found'])
    pause(2)
end

if isunix
    
    misc.username           = 'alejandro';
    misc.Databases          = ['~/Documenten/Databases' delim];
    misc.db_speechmaterials_local = '/media/Elements/orl-wrk-0089/Documenten/Speech_material_from_x-drive/';
    misc.db_speechmaterials = ['/home/' misc.username '/x/speechmaterials/'];
    misc.ex_APEX_results    = ['~/Documenten/Meas/Meas/Experiments/Results_XML/'];
    misc.lx_Templates       = ['~/Documenten/Documenten-TUe/01-Text/00-Templates' delim];
    misc.lx_Text            = ['~/Documenten/Documenten-TUe/01-Text/05-Doc-TUe'   delim];
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
    misc.Databases          = ['D:\Databases' delim];
    misc.SVN_KUL            = ['D:\SVN-KU-Leuven\alejandro' delim];
    misc.db_speechmaterials_local = ['E:\orl-wrk-0089\Documenten\Speech_material_from_x-drive']; % needs my external hard drive connected (E-drive)
    misc.db_speechmaterials = ['D:\Databases\dir03-Speech'];
    misc.db_fastl2007       = ['D:\Documenten-TUe\10-Referenties\02-Mijn-boeken\Fastl2007-psychoacoustics' delim 'Sound' delim];
    misc.ex_APEX_results    = [misc.SVN_KUL 'Meas' delim 'Experiments' delim 'Results_XML' delim]; % work at KUL
    misc.lx_Presentations   = ['D:\Documenten-TUe\01-Text\70-Presentaties-TUe' delim];
    misc.lx_Templates       = ['D:\Documenten-TUe\01-Text\00-Templates' delim];
    misc.lx_Text            = ['D:\Documenten-TUe\01-Text\05-Doc-TUe'   delim];
    misc.outputs            = ['D:\MATLAB\Output' delim];
    misc.praat              = ['C:\praat5376_win32\praatcon.exe'];
    
end

misc.db_audacity        = ['C:\Users\aosses\Documents\db_audacity' delim];

misc.db_voice_of_dragon = [misc.Databases 'dir01-Instruments' delim 'Voice-of-dragon' delim];
misc.db_HRIR_Oldenburg  = [misc.Databases 'dir02-HRTFs'       delim 'Oldenburg'       delim];
misc.db_Fastl2007       = [misc.Databases 'dir04-Psychoacoustics' delim 'Fastl-and-Zwicker-2007' delim '03-Extracted-files' delim];
misc.db_calfiles        = [misc.Databases 'dir04-Psychoacoustics' delim 'Fastl-and-Zwicker-2007' delim '03-Extracted-files' delim]; % TMP dir
misc.language           = 'NL'; % other possibilities: 'EN'
     
if (nargin==1)
    if (~isfield(misc,type))
        error('Invalid type');
    end
    misc=misc.(type);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

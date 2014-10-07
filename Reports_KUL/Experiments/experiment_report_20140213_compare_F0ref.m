function experiment_report_20140213_compare_F0ref(info)
% function experiment_report_20140213_compare_F0ref(info)
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localpaths = Get_paths([],1); % addpath('~/Documenten/MATLAB/MATLAB_svn/Utility/')

bPart = [1 0 0];

info.bAssess    = 0;
info.bAnalyse   = 1;
info.isPaper    = 1;
info.bSave      = 0;

if info.bSave
    disp([mfilename '.m: Figures and other files will be stored and might replace other files. Press any button to continue'])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause()
else
    disp([mfilename '.m: Figures and other files won''t be stored...change bSave to 1 if want to save results to file...'])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause(3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.F0reference        = 'laryngeal';
info.F0max              = 400;

info.bCleanSpeech       = 1; % if 0: CP810 simulations will be used
info.bPlot              = 0;
info.speaker            = 'sb'; % Female 
info.results_dir_name   = 'results';
info.figures_dir_name   = 'figures';



if bPart(1) == 1 & info.isPaper
    
    info = experiment_report_20140228_Physical_validation(info);
    
end

info.speaker            = 'rl'; % Male
if bPart(1) == 1
    experiment_report_20140228_Physical_validation(info);
end

% Analysis as in Vandali 2011 for female and male speakers
error_tot_voiced_ref = experiment_report_20140228_F0_as_Vandali(info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
info.F0reference        = 'praat';
info.F0max              = 400;
  
info.bCleanSpeech       = 1; % 0 = Then it will use CP810 simulations
info.speaker            = 'sb'; % Female 
info.results_dir_name   = 'results-praat-clean';
info.figures_dir_name   = 'figures-praat-clean';

if bPart(2) == 1
    info.bSave = 1;
    info = experiment_report_20140228_Physical_validation(info);
end

info.speaker            = 'rl';

if bPart(2) == 1
    info.bSave = 1;
    info = experiment_report_20140228_Physical_validation(info);
end

% % Analysis as in Vandali 2011 for female and male speakers
% error_tot_voiced_ref = experiment_report_20140228_F0_as_Vandali(info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST-f Evaluation

info.F0reference        = 'praat';
info.F0max              = 400;

info.bCleanSpeech       = 0;
info.speaker            = 'sb'; % Female 
info.results_dir_name   = 'results-praat-CP810';
info.figures_dir_name   = 'figures-praat-CP810';

info_LIST = info;
info_LIST.bAssess = 0;
info_LIST.bAnalyse = 1;
info_LIST.F0max = 300;
info_LIST = experiment_report_20140228_Physical_validation_LIST(info_LIST);

error_tot_voiced_ref = experiment_report_20140228_F0_as_Vandali(info_LIST);

h = [];

h(end+1) = Plot_errorrate_F0mod(error_tot_voiced_ref);
filename = [info_LIST.figures 'ac-error-rate-LIST'];
if info_LIST.bSave
    Saveas(h(end), filename);
end

% End of LIST-f evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bPart(3) == 1
    info.bSave = 1;
    info = experiment_report_20140228_Physical_validation(info);
end

info.speaker            = 'rl';

if bPart(3) == 1
    info.bSave = 1;
    info = experiment_report_20140228_Physical_validation(info);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis as in Vandali 2011 for female and male speakers
error_tot_voiced_CP810 = experiment_report_20140228_F0_as_Vandali(info);

error_tot = [error_tot_voiced_ref; error_tot_voiced_CP810];

% This is going to be for my thesis
h = Plot_errorrate_F0mod_vs_eTone(error_tot);

info.figures = '~/Documenten/fda_eval/figures/';

if info.bSave
    filename = [info.figures 'ac-vs-eTone2011'];
    Saveas(h, filename);
end

h(end+1) = Plot_errorrate_F0mod(error_tot_voiced_ref);
filename = [info.figures 'ac-error-rate'];
if info.bSave
    Saveas(h(end), filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
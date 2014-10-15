function paper_IJA_2014
% function paper_IJA_2014
%
% Created in    : 2014
% Last update on: 25/08/2014
% Last use on   : 25/08/2014
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Commit:
%   experiment_report_20130131_VlMatrix
%   Check output dir fda evaluation

fprintf('Make sure analysis of PB and LIST-f databases are in the proper directory...\n');
pause(2)

bDoErr      = 1; % Figure 2
bDoPR       = 0;
bDoSpeech   = 0;
options.bSave = 1;

% If bDoPR
if isunix
    infodir.export_folder   = '~/Documenten/LaTeX_Docs/paper/data_new/'; % local folder Alejandro
    info.figures_folder  = '~/Documenten/LaTeX_Docs/paper/figures_new/'; % local folder Alejandro
else
    infodir.export_folder   = [Get_TUe_paths('outputs') 'IJA_data' delim]; % local folder Alejandro
    info.figures_folder  = [Get_TUe_paths('outputs') 'IJA_figures' delim]; % local folder Alejandro
end

try
    % cd('~/repo/alejandro/documenten/paper/'); % if Tom's computer
    cd('~/Documenten/LaTeX_Docs/paper/')
    dirtree.paper_folder    = cd;
catch
    dirtree.paper_folder    = [Get_TUe_paths('outputs') 'paper' delim]; % uigetdir('~/','Select an output directory for physical validation');
    Mkdir(dirtree.paper_folder);
end

% bDoSpeech
options.dest_folder = [dirtree.paper_folder 'IJA_figures_speech' delim]; 

% bDoPR and bDoErr
info.export_folder   = [dirtree.paper_folder 'data_new' delim]; % local folder Alejandro
Mkdir(info.figures_folder); % Folder for F0 error plots
info.bSave           = options.bSave;

Mkdir(options.dest_folder);
Mkdir(info.export_folder);



%% Figure 2, 3
if bDoErr

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % F0 performance
    experiment_report_20140402_compare_F0ref(info); % This script replaces % plottingErrors; % (20 sentences analysis)
    % options.bPlot = 1;
    % experiment_report_20140228_F0_as_Vandali(options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

if isunix
    addpath([info.svn.MATLAB 'Simulink_XPC/']); % Running_XPC_off_post2.m
    close all; 
end

try
    % cd('~/repo/alejandro/documenten/paper/'); % if Tom's computer
    cd('~/Documenten/LaTeX_Docs/paper/')
    dirtree.paper_folder    = cd;
catch
    if isunix
        dirtree.paper_folder    = uigetdir('~/','Select an output directory for physical validation');
    else
        dirtree.paper_folder = Get_TUe_paths('outputs');
    end
end

dirtree.root_location   = [dirtree.paper_folder delim 'validation' delim];

%%
if bDoPR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PR Results (Figure 8):
    infodir.bSave           = options.bSave;
    f0mod_pr_results(infodir); % Programmed by Tom Francart
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
if bDoSpeech
    options.dest_folder = [dirtree.paper_folder delim 'figures' delim];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bDoSpeechTests = 1;
    % bDoLT       = 1;
    % bDoLL       = 1;
    % bDoMT       = 1;
    % bDoStats    = 0; % Not tested
    %
    % Tested on 25/08/2014: (Save not tested)
    experiment_report_20131127_pooled(options); % Lilliput, LIST, VlMatrix results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    experiment_report_20130131_VlMatrix(options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
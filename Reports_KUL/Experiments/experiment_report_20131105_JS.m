function experiment_report_20131105_JS
% function experiment_report_20131105_JS
%
% Analysing Julia Schoolmeesters' data, Lilliput, VMatrix
%
% Dependencies:
%   figLTScores_xPC
% 	quick_LL_one_subject
%   quick_MT_one_subject
%
%   Individual results are extracted from global txt-files:
%
% Results for LIST      = 'ci-Julia_Schoolmeesters/CI-Julia-LT-adaptive.txt'
% Results for VMatrix   = 'CI-Pooled-MT-fixedSNR.txt'
% Results for Lilliput  = 'CI-Pooled-Word-LL.txt'
%
% Programmed by Alejandro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Statistics/') % Path for Get_mean
 
hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/';
subject_folder  =  'ci-Julia_Schoolmeesters/';
dest_folder     = [hoofd_folder 'Experiment_reports/20131105-ci-Julia_Schoolmeesters/Figures/'];
% dest_folder     = [hoofd_folder subject_folder 'samenvatting-20131105/Figures/'];

directories2check = {[subject_folder '20131105-LT/']};

% SubjectName     = 'CIs';
subjects        = 17; % Julia Schoolmeesters

file_LIST_adaptive  = [hoofd_folder 'Results_XML/' subject_folder '/CI-Julia-LT-adaptive.txt'];
file_MT_fixedSNR    = [hoofd_folder 'Results_XML/' 'CI-Pooled-MT-fixedSNR.txt']
file_LL_fixedSNR    = [hoofd_folder 'Results_XML/' 'CI-Pooled-Word-LL.txt'];
 
bSave = 1;
bPlot = 1;
bDoSpeechTests = 1;
nSessie = 2; 
 
if bDoSpeechTests == 1
    
    stPlot.xTick            = 1;
    stPlot.bPlotIndividual  = 1;
    stPlot.YGrid            = 'on';
    stPlot.yLim             = [-2 6];
    [xx pLT]                = figLTScores_xPC(file_LIST_adaptive, bPlot, stPlot);

    
    if bSave
        quick_staircases(directories2check, dest_folder, hoofd_folder);
    end
    
    
    if bSave
        saveas(gcf,[dest_folder,'LT-adaptive-results'],'epsc');
    end
    
    [h p]       = figLLScores_xPC(file_LL_fixedSNR);
    stHandle = quick_LL_one_subject(p, subjects, nSessie);
 
    if bSave
        saveas(stHandle.h1     ,[dest_folder,'Lilliput-one-subject-report'   ,num2str(nSessie)],'epsc');
        % saveas(stHandle.hPooled,[dest_folder,'Lilliput-one-subject-up2report',num2str(nSessie)],'epsc'); % pooled data = individual data
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATRIX: 

    stPlot.bPlotIndividual = 1;
    [handleFig pMT] = figMTScores_xPC(file_MT_fixedSNR,stPlot);

    close;
    close; % close the 2 generated plots

    % Matrix word-scores (plot for 1 individual) 
 
    stHandle = quick_MT_one_subject(pMT, subjects, nSessie);

    if bSave
        saveas(stHandle.h1     ,[dest_folder,'Vlaamse-matrix-one-subject-report'   ,num2str(nSessie)],'epsc');
        % saveas(stHandle.hPooled,[dest_folder,'Vlaamse-matrix-one-subject-up2report',num2str(nSessie)],'epsc'); % pooled data = individual data
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
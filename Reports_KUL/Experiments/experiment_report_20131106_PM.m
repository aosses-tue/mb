function experiment_report_20131106_PM
% function experiment_report_20131106_PM
%
% Analysing Patrick Meul's data, Lilliput, VMatrix, PR (Loudness Balancing 
% as done in previous session)
%
% Dependencies:
%   figLTScores_xPC
%   figLLScores_xPC
%
%   Individual results are extracted from global txt-files:
%
% Results for VMatrix   = 'CI-Pooled-MT-fixedSNR.txt'
% Results for Lilliput  = 'CI-Pooled-Word-LL.txt'
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Statistics/') % Path for Get_mean

hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
dest_folder     = [hoofd_folder 'ci-Patrick_Meul/samenvatting-20131106/Figures/'];

SubjectName     = 'CIs';
subjects        = 13; % Patrick Meul

file_MT_fixedSNR    = [hoofd_folder 'CI-Pooled-MT-fixedSNR.txt']
file_LL_fixedSNR    = [hoofd_folder 'CI-Pooled-Word-LL.txt'];

% bPlot = 1;
bSave = 1;
bDoSpeechTests = 1;
% bDoPR = 1;
nSessie = 2; 

if bDoSpeechTests == 1
     
    [h p]       = figLLScores_xPC(file_LL_fixedSNR);
    stHandle = quick_LL_one_subject(p, subjects, nSessie);
 
    if bSave
        saveas(stHandle.h1     ,[dest_folder,'Lilliput-one-subject-report'   ,num2str(nSessie)],'epsc');
        saveas(stHandle.hPooled,[dest_folder,'Lilliput-one-subject-up2report',num2str(nSessie)],'epsc');
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
        saveas(stHandle.hPooled,[dest_folder,'Vlaamse-matrix-one-subject-up2report',num2str(nSessie)],'epsc');
    end

end

stPlot = [];
stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE'};
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pitch ranking

stSubject.truncate_data = 0;

stSubject.firstname = 'ci-Patrick';
disp('When asked, press 1')
disp('  ')
disp('  ')

[aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);

if bSave
    saveas(gcf,[dest_folder,'PR-results-PM-old'],'epsc');
end

disp('When asked, press 2')
disp('  ')
disp('  ')

[aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);

if bSave
    saveas(gcf,[dest_folder,'PR-results-PM-retest'],'epsc');
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
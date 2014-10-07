function experiment_report_20131127_WD
% function experiment_report_20131127_WD
%
% Analysing Wouter David's data, Lilliput, VMatrix, PR, LB
%
% Dependencies:
%   figLTScores_xPC
%   figLLScores_xPC
%
%   Individual results are extracted from global txt-files:
%
% Programmed by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Statistics/') % Path for Get_mean

info = getpaths;

SubjectName     = 'CIs';
subjects        = 15; % Wouter David

session_date    = '20131127';
hoofd_folder    = [info.svn.Meas 'Experiments/'];
test_folder     = [hoofd_folder 'Results_XML/'];
result_folder   = info.result_folder;
subject_folder  =  'ci-Wouter_David/';

stSubject       = Get_subject_info(subject_folder);
% subjects        = stSubject.nSubject;

dest_folder     = [info.experiment_report session_date '-' subject_folder 'Figures/'];

% bPlot = 1;
bSave   = 1;
bDoSpeechTests = 0;
bDoPR   = 1;

if bDoSpeechTests == 1
     
    %[h p]       = figLLScores_xPC(file_LL_fixedSNR);
    %stHandle = quick_LL_one_subject(p, subjects, nSessie);
 
    if bSave
        %saveas(stHandle.h1     ,[dest_folder,'Lilliput-one-subject-report'   ,num2str(nSessie)],'epsc');
        %saveas(stHandle.hPooled,[dest_folder,'Lilliput-one-subject-up2report',num2str(nSessie)],'epsc');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATRIX: 

    %stPlot.bPlotIndividual = 1;
    %[handleFig pMT] = figMTScores_xPC(file_MT_fixedSNR,stPlot);

    %close;
    %close; % close the 2 generated plots

    % Matrix word-scores (plot for 1 individual) 
 
    %stHandle = quick_MT_one_subject(pMT, subjects, nSessie);

    if bSave
        % saveas(stHandle.h1     ,[dest_folder,'Vlaamse-matrix-one-subject-report'   ,num2str(nSessie)],'epsc');
        % % saveas(stHandle.hPooled,[dest_folder,'Vlaamse-matrix-one-subject-up2report',num2str(nSessie)],'epsc');
    end

end


if bDoPR
    stPlot = [];
    stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pitch ranking

    stSubject.truncate_data = 0;

    stSubject.firstname = 'ci-Wouter';
    disp('When asked, press 1')
    disp('  ')
    disp('  ')

    stPlot.scalerate = 0.45;
    
    [aceData, f0mData, p] = figPRScores_xPC(test_folder,{'UW'}, stPlot, stSubject);
    %%%%%%%%%%%%%%%%%%%%%
    % Excluding non-measured data:
    % UW_262
    aceData(17:20) = NaN;
    f0mData(17:20) = NaN;
    
    % UW_131
    f0mData(5:8) = NaN;
    %%%%%%%%%%%%%%%%%%%%%
    figPRScores_xPC_only_plot(aceData,f0mData,stPlot,p,stSubject);
    
    if bSave
        saveas(gcf,[dest_folder,'PR-results-WD'],'epsc');
    end
    
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

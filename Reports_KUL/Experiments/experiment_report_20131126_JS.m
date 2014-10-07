function experiment_report_20131126_JS
% function experiment_report_20131126_JS
%
% Analysing Julia Schoolmeesters' data, Lilliput, VMatrix, PR, LB
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
% Programmed by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Statistics/') % Path for Get_mean

info = getpaths;

session_date    = '20131126';
hoofd_folder    = [info.svn.Meas 'Experiments/'];
result_folder   = info.result_folder;
subject_folder  =  'ci-Julia_Schoolmeesters/';

stSubject       = Get_subject_info(subject_folder);
subjects        = stSubject.nSubject; % Julia Schoolmeesters

dest_folder     = [info.experiment_report session_date '-' subject_folder 'Figures/'];

file_LISTf      = [result_folder delim subject_folder delim session_date '-LISTf' delim session_date '-CI-JS-LISTf'];
file_MT         = [result_folder delim subject_folder delim session_date '-MT'    delim session_date '-CI-JS-VlMatrix.txt'];
file_LL         = [result_folder delim subject_folder delim                             session_date '-CI-JS-LL.txt'];

bSave   = 0;
bDoSpeechTests = 1;
bDoLIST = 1;
bDoMT   = 0;
bDoLL   = 0;

if bDoSpeechTests == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LIST-f
    
    if bDoLIST
        
        if bSave
            quick_staircases(directories2check, dest_folder, hoofd_folder);
        end

        stPlot.TitleHead = 'LIST-f';
        stPlot.YLabel = 'SRT (dB)';
        stPlot.XLabel = {'before'};
        stPlot.SeriesLabel = {'xPC F0mod','xPC ACE'};
        stPlot.bPlotIndividual = 1;
        stPlot.YGrid            = 'on';
        stPlot.yLim             = [-7 7];
        stPlot.xLim             = [0 3];

        [xx pLTf_b] = figLTScores_xPC([file_LISTf '-before.txt'], 1, stPlot);

        if bSave == 1
            saveas(gcf,[dest_folder, 'LISTfbefore.eps'],'epsc');
        end

        stPlot.XLabel = {'after'};

        [xx pLTf_a] = figLTScores_xPC([file_LISTf '-after.txt' ], 1, stPlot);   

        if bSave == 1
            saveas(gcf,[dest_folder, 'LISTfafter.eps'],'epsc');
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lilliput, just subject S17
    
    if bDoLL
        
        [h p]           = figLLScores_xPC(file_LL);
                
        for i = 1:length(subjects) % Individual pooled results. One plot per subject
            stHandle = quick_LL_one_subject(p, subjects(i));

            if bSave
                saveas(stHandle.h1     ,[dest_folder,'Lilliput-one-subject-report-S'   ,num2str(subjects(i))],'epsc');
                close
            end
        end
         
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if bDoMT
        
        directory_VlMatrix = [hoofd_folder 'Results_XML/' subject_folder delim session_date '-MT' delim];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % training for xPC F0mod:
        % F0mod_results = dir([directory_VlMatrix 'dubbel*F0m-*' stSubject.Initials '.apr']);
        F0mod_results = dir([directory_VlMatrix '*lijst_*F0m-*' stSubject.Initials '.apr']);

        resF0m.SNR = [];
        for i = 1:length(F0mod_results)
            [xx xx info] = parseapex_VlMatrix([directory_VlMatrix F0mod_results(i).name]);

            if str2num(info.SNR) > 50
                info.SNR = 99;  
            else
                info.SNR = str2num(info.SNR);
            end

            resF0m.wscores(i)= info.score_words;        
            resF0m.wscores1(i) = info.reliability(1);
            resF0m.wscores2(i) = info.reliability(2);
            resF0m.SNR(i)    = info.SNR;
        end
        
        idx = find(resF0m.SNR == 10);
        
        figure
        plot(1:length(idx), resF0m.wscores(idx), 'bo','LineWidth',4), hold on, grid on
        Title('VlMatrix with SNR = 10 dB, xPC F0mod')
        Legend('Scores F0mod')
        Xlabel('Order presented')
        Ylabel('Word score [%]')
        ylim([20 100])
        xlim([0 length(idx)+1])
             
        if bSave == 1
            saveas(gcf,[dest_folder, 'LearningMT.eps'],'epsc');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % xPC ACE: normally we will not need to analyse ACE learning
        ACE_results = dir([directory_VlMatrix '*ACE-*' stSubject.Initials '.apr']);

        resACE.SNR = [];
        for i = 1:length(ACE_results)
            [xx xx info] = parseapex_VlMatrix([directory_VlMatrix ACE_results(i).name]);

            if str2num(info.SNR) > 50
                info.SNR = 99;  
            else
                info.SNR = str2num(info.SNR);
            end

            resACE.wscores(i)= info.score_words;        
            resACE.wscores1(i) = info.reliability(1);
            resACE.wscores2(i) = info.reliability(2);
            resACE.SNR(i)    = info.SNR;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [xxx pMT] = figMTScores_xPC(file_MT);
        handle = quick_MT_one_subject(pMT, subjects);

        if bSave == 1
            saveas(handle.h1,[dest_folder, 'MT-one-subject-report-' stSubject.Initials '.eps'],'epsc');
        end
        
    end
   
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
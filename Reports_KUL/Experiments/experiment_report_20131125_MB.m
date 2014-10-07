function experiment_report_20131125_MB
% function experiment_report_20131125_MB
%
% Analysing Maria Brughmans' data: LIST-f, VlMatrix, PR (only UW 131 Hz)
%
% Dependencies:
%   figLTScores_xPC
%   figLLScores_xPC
%
%   Individual results are extracted from global txt-files:
%
% Results for VMatrix   = 'CI-Pooled-MT-fixedSNR.txt'
%
% Programmed by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Statistics/') % Path for Get_mean
info = getpaths;
List_of_files   = {     'setupPlotConf.m'};
Check_dependencies_ExpORL(List_of_files, info.username);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session_date    = '20131125';
hoofd_folder    = [info.svn.Meas 'Experiments/'];
result_folder   = info.result_folder;
subject_folder  =  'ci-Maria_Brughmans/';

stSubject       = Get_subject_info(subject_folder);
subjects        = stSubject.nSubject; % Maria Brughmans

dest_folder     = [info.experiment_report session_date '-' subject_folder 'Figures/'];

file_LISTf      = [result_folder delim subject_folder delim session_date '-LISTf' delim session_date '-CI-MB-LISTf'];
file_MT         = [result_folder delim subject_folder delim session_date '-MT'    delim session_date '-CI-MB-VlMatrix.txt'];

% For plotting staircases:
directories2check   ={['Results_XML' delim subject_folder delim session_date '-LISTf' delim]};

bSave   = 0;
bDoSpeechTests = 1;
bDoLIST = 1;
bDoMT   = 0;
bDoPR   = 0;

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
        stPlot.SeriesLabel = {'xPC F0mod',''};
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
    
    if bDoMT
        
        directory_VlMatrix = [hoofd_folder 'Results_XML/' subject_folder delim session_date '-MT' delim];
        
        % xPC F0mod:
        F0mod_results = dir([directory_VlMatrix 'dubbel*F0m-*' stSubject.Initials '.apr']);

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
        
        % xPC ACE: % Only 2 basislisjten
        ACE_results = dir([directory_VlMatrix 'basis*ACE-*' stSubject.Initials '.apr']);

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
        [xxx pMT] = figMTScores_xPC(file_MT);
        handle = quick_MT_one_subject(pMT, subjects);

        if bSave == 1
            saveas(handle.h1,[dest_folder, 'MT-one-subject-report-' stSubject.Initials '.eps'],'epsc');
        end
        
    end
   
end

if bDoPR
    
    stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE'};
    disp('   ')
    disp([mfilename '.m: type 2 when requested'])
    disp('   ')
    
    % Height of screen > 1050 pixels
    stPlot.scalerate = 0.5; % Comment this if the screen height is greater than 1050 pixels
    [aceData, f0mData]      = figPRScores_xPC(result_folder,{'UW'}, stPlot, stSubject);
    if bSave
        saveas(gcf,[dest_folder,'PR-results-MB-old'],'epsc');
    end
    
    disp('   ')
    disp([mfilename '.m: type 3 when requested'])
    disp('   ')
    
    [aceData, f0mData]      = figPRScores_xPC(result_folder,{'UW'}, stPlot, stSubject);
    if bSave
        saveas(gcf,[dest_folder,'PR-results-MB-new'],'epsc');
    end
    
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
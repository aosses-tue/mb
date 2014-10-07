function experiment_report_20131122_JL
% function experiment_report_20131122_JL
%
% Analysing Jan Leys' data, LIST, PR. LB (re-used)
%
% Dependencies:
%       
%   KeepCapitalLetters
%   figLTScores_xPC
%   quick_staircases
%
%   Results for LIST adaptive = 'ci-Jan_Leys/CI-Jan-LT-adaptive.txt'
%   Results for LIST-f adaptive (training) = 'ci-Jan_Leys/20131128-CI-Jan-LISTf.txt'
%   Results for LIST fixed = handwritten
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Statistics/') % Path for Get_mean
 
hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/';
test_folder     = [hoofd_folder 'Results_XML/'];
subject_folder  =  'ci-Jan_Leys';

stSubject.Initials = KeepCapitalLetters(subject_folder); % JL
stSubject.firstname = 'ci-Jan';
stSubject.truncate_data = 0;

file_LISTf      = [hoofd_folder 'Results_XML/' subject_folder '/20131128-LIST-f/20131128-CI-Jan-LISTf'];
file_LISTm      = [hoofd_folder 'Results_XML/' subject_folder '/20131128-LIST-m/20131128-CI-Jan-LISTm'];

% For plotting staircases:
directories2check   ={['Results_XML' delim subject_folder delim '20131122-LIST-f' delim], ...
                      ['Results_XML' delim subject_folder delim '20131122-LIST-m' delim]};

% output path for figures:
dest_folder         = [hoofd_folder 'Experiment_reports/20131122-' subject_folder '-training/Figures/'];
directory_VlMatrix  = [hoofd_folder 'Results_XML/' subject_folder '/20131122-MT/'];

% SubjectName     = 'CIs';
% subjects        = 14; % Jan Leys
 
bPlot   = 1;
bSave   = 0;
bDoSpeechTests = 1;

if bDoSpeechTests == 1
     
    if bSave
        quick_staircases(directories2check, dest_folder, hoofd_folder);
    end
    
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

    idx = find(resF0m.SNR == 5);
    
    plot(1:length(idx), resF0m.wscores(idx), 'bo','LineWidth',4), hold on, grid on
    Title('VlMatrix with SNR = 5 dB, xPC F0mod')
    Legend('Scores F0mod')
    Xlabel('Order presented')
    Ylabel('Word score [%]')
    ylim([20 100])
    xlim([0 length(idx)+1])
    
    if bSave == 1
        saveas(gcf,[dest_folder, 'LearningMT.eps'],'epsc');
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
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    stPlot.TitleHead = 'LIST-m';
    [xx pLTm_b] = figLTScores_xPC([file_LISTm '-before.txt'], 1, stPlot);
    
    if bSave == 1
        saveas(gcf,[dest_folder, 'LISTmbefore.eps'],'epsc');
    end
    
    stPlot.XLabel = {'after'};
    
    [xx pLTm_a] = figLTScores_xPC([file_LISTm '-after.txt' ], 1, stPlot);   
    
    if bSave == 1
        saveas(gcf,[dest_folder, 'LISTmafter.eps'],'epsc');
    end
    
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
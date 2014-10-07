function experiment_report_20131108_JL
% function experiment_report_20131108_JL
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
%   Results for LIST fixed = handwritten
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Statistics/') % Path for Get_mean

% hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/';
hoofd_folder    = '/home/tom/temp/alejandro/Meas/Experiments/';
test_folder     = [hoofd_folder 'Results_XML/'];
subject_folder  =  'ci-Jan_Leys/';

stSubject.Initials = KeepCapitalLetters(subject_folder); % JL
stSubject.firstname = 'ci-Jan';
stSubject.truncate_data = 0;

file_LIST_adaptive  = [hoofd_folder 'Results_XML/' subject_folder '/CI-Jan-LT-adaptive.txt'];
directories2check   ={['Results_XML/' subject_folder '20131003-LT/']};

dest_folder     = [hoofd_folder 'Experiment_reports/20131108-' subject_folder 'Figures/'];

% SubjectName     = 'CIs';
% subjects        = 14; % Jan Leys

bPlot   = 1;
bSave   = 0;
bDoSpeechTests = 1;
bDoPR   = 0;
FontSize = 14;

if bDoSpeechTests == 1

    
    stPlot.xTick            = 1;
    stPlot.bPlotIndividual  = 1;
    stPlot.YGrid            = 'on';
    stPlot.yLim             = [-6 6];
    [xx pLT]                = figLTScores_xPC(file_LIST_adaptive, bPlot, stPlot);
    
    if bSave
        saveas(gcf,[dest_folder,'LT-adaptive-results'],'epsc');
    end
    
    if bSave
        quick_staircases(directories2check, dest_folder, hoofd_folder);
    end
    
    % Automate the following:
    PsychoACE = [   0.4 0.6 0.9 0.8; ... % scores (0-1)
                    0   2   3   5];      % SNR
                
	PsychoF0m = [   0.1 0.3 0.3 0.5 0.5 0.6 0.6; ...
                    0   2   3   4   5   6   7];
    Xmin = -2;
    Xmax = 8;
    
    figure,
	plot(PsychoACE(2,:), PsychoACE(1,:)*100, 'bo', 'LineWidth',5), grid on, hold on
	ylim([0 100])
    xlim([Xmin 8])
    
    plot(PsychoF0m(2,:), PsychoF0m(1,:)*100, 'ro', 'LineWidth',5)
    hLabel = ylabel('Sentence scores [%]'); set(hLabel,'FontSize',FontSize)
    hLabel = xlabel('SNR (dB)');        set(hLabel,'FontSize',FontSize)
    hLabel = legend('xPC ACE','xPC F0mod');     set(hLabel,'FontSize',FontSize)
    
    hLabel = gca; 
    set(hLabel,'FontSize',FontSize)
    set(hLabel,'XTick',[Xmin:Xmax])
    
    plot(PsychoACE(2,:), PsychoACE(1,:)*100, 'b','LineWidth',0.5)
    plot(PsychoF0m(2,:), PsychoF0m(1,:)*100, 'r','LineWidth',0.5)
    
    plot([-8 8],[50 50],'k--','LineWidth',2)
    
    if bSave
        saveas(gcf,[dest_folder,'LT-fixed-results'],'epsc');
    end
    
end

if bDoPR
    stPlot = [];
    stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pitch ranking

    disp('When asked, press 1')
    disp('  ')
    disp('  ')

    [aceData, f0mData]      = figPRScores_xPC(test_folder,{'UW'}, stPlot, stSubject);

    if bSave
        saveas(gcf,[dest_folder,'PR-results-' stSubject.Initials '-old'],'epsc');
    end

    disp('When asked, press 2')
    disp('  ')
    disp('  ')

    [aceData, f0mData]      = figPRScores_xPC(test_folder,{'UW'}, stPlot, stSubject);

    if bSave
        saveas(gcf,[dest_folder,'PR-results-' stSubject.Initials '-retest'],'epsc');
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
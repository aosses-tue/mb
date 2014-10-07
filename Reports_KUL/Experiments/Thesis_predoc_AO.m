function Thesis_predoc_AO(dest_folder)
% function Thesis_predoc_AO(dest_folder)
%
% Generating plots used in my predoc thesis
%
% dest_folder - folder where the figures will be saved
%
% Tested on Windows 7 just analysing Music tests (set the rest of boolean
% variables to 0)
% 
% % Example, Wekelijkse update 30:
%   dest_folder = '/home/alejandro/Documenten/LaTeX_Docs/predoc_plan_wekelijkse_update/wu031_2013_09_03/Figures'
%   Thesis_predoc_AO(dest_folder)
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven 2014
% Created on:
% Last updated on: 12/06/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isunix
    
    List_of_files = {   'single_channel_sequence.m', ...
                        'thesis_chapter1fig05_EHPitchMech.m'};

    [AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDoSpeechTests  = 0;
bDoMusicTests   = 1;
bPrepareSPSS    = 0;
bDoOtherPlots   = 0;
bSave           = 0;

disp([mfilename '.m: it is possible to customise which experiments are going to be processed... (edit this file if desired)'])
pause(1)

if nargin == 0
    if isunix
        dest_folder = '/home/alejandro/Documenten/LaTeX_Docs/predoc_plan/proefschrift/Figures-auto_new/';
        result_folder = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
    else
        dest_folder = Get_TUe_paths('outputs');
        result_folder = Get_TUe_paths('ex_APEX_results');
    end
end

h = [];

if bDoOtherPlots
    h(end+1) = thesis_chapter1fig05_EHPitchMech;
    if bSave
        saveas(h(end),[dest_folder,'EH-Pitch'],'epsc')
    end
end

if bDoSpeechTests
    h(end+1) = figLLScores_xPC;
    if bSave
        saveas(h(end),[dest_folder,'Word-results'],'epsc')
    end

    [h(end+1) pLT] = figLTScores_xPC;
    if bSave    
        saveas(h(end),[dest_folder,'Sentence-results'],'epsc')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This test was not included in the thesis
    [h(end+1), pLT_old] = figLTScores_xPC('NH-Sentence-recognition-test-4-subjects.txt');
    [h(end+1), pLT_new] = figLTScores_xPC('NH-Sentence-recognition-retest2.txt'); % retest-4subjects
    stPlot.TitleHead = ['LIST material (FIR EE)'];
    stPlot.xTickLabel       = {'old','new'};
    h(end+1) = PlotMeans_default(stPlot,pLT_old.aceData, pLT_old.f0mData, pLT_new.aceData, pLT_new.f0FIREEData );
    saveas(h(end),[dest_folder,'Sentence-results+retest'],'epsc')
end

if bDoMusicTests
    
    [aceData, f0mData] = figPRScores_xPC(result_folder,{'UW'});
    % Then type (when requested) the following:
    %   2
    %   2
    %   [1 1 1 1 1 1 1] 
    h(end+1) = gcf;
    if bSave
        saveas(h(end)  ,[dest_folder,'PR'],'epsc')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preparing SPSS:
    if bPrepareSPSS
        convertPR2SPSS(aceData,f0mData);
    end
    h(end+1:end+2) = figMCScores_xPC(result_folder,{'UW'});
    % Then type (when requested) the following:
    %   3
    %   2
    %   [1 0 1 1 1 1 1] % To exclude Anneke's data (incomplete by 08-08-2013)
    if bSave
        saveas(h(end-1),[dest_folder,'MC-per-interval'],'epsc')
        saveas(h(end)  ,[dest_folder,'MC-per-sequence'],'epsc')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
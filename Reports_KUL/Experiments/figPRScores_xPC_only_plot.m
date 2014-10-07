function figPRScores_xPC_only_plot(aceData, f0mData, stPlot, p, stSubject)
% function figPRScores_xPC_only_plot(aceData, f0mData, stPlot, p, stSubject)
%
% Adapted from: thesis_chapter3PRScores
%
% Figure 3.6, Page 58, Matthias' Thesis
%
% stPlot    - structure containing Plotting options using Matthias' GUI        
% stSubject - structure with Subject's info, some fields:
%       bRoving: when set to 1 plots the roving values across registers
%
% Dependencies:
%   readPRData.m
%
% Examples:
%       close all, 
%       clc, 
%       figPRScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'CLARINET'})
%
%       close all, 
%       clc, 
%       [aceData, f0mData] = figPRScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'UW'});
%
% Programmed by Alejandro, adapted from Matthias' script thesis_chapter3PRScores
%
% Things to Implement:
%       1. Automatic search of Experiments with ACE and F0mod strategies for 
%          all the found subjects (implement this on this script and figMCScores_xPC
%       2. To check significance level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentDirectory = cd;
List_of_files = {   'get_tones_for_pitch_ranking.m', ... %'confidence_interval.m', ... %                     'get_significance_level.m', ...
                    'get_significance_level.m', ...
                    'eval_results.m', ...
                    'note2freq.m', ...
                    'plot_PR_results.m', ...
                    'setupPlotConf.m'};
                
Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    instr{1} = '';
end

if nargin < 4
    stSubject = [];
end

stPlot      = Ensure_field(stPlot,'scalerate',1);
stSubject   = Ensure_field(stSubject,'bRoving',0);

regs        = {'G#2', 'C3', 'E3', 'G#3', 'C4'};

nRegisters  = length(regs);

disp(['Running script: ' mfilename])
disp('Please check that all the experiments for the registers are in your computer')
disp(['If don''t this script will fail. You are trying to load: ' num2str(nRegisters)])

% [aceData, f0mData, p, stRove] = readPRData(dirExpRes, regs, instr{1}, stSubject);
[nSubjects bb]      = size(aceData);
disp(['Number of subjects taken into account: ', num2str(nSubjects)])

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nRepetitions        = p.nRepetitions; 
nStrategies         = p.nStrategies;
nSemitones          = p.nSemitones;

stPlot.nPlots       = p.nReferences;
stPlot.nCols        = 2;
stPlot.nRows        = ceil(p.nReferences/stPlot.nCols);

stPlot.xGap         = 15*stPlot.scalerate;
stPlot.yGap         = 80*stPlot.scalerate;
stPlot.fntsz        = 16*stPlot.scalerate;
stPlot.markerSize   = 12*stPlot.scalerate;
stPlot.margins      = [70 70 50 90]*stPlot.scalerate;
% figPos              = [0 0 1024 225*nRows]; stPlot.figPos       = figPos; % 225 each row
stPlot              = Ensure_field(stPlot, 'figPos',[0 0 1024 350*stPlot.nRows]*stPlot.scalerate); % 225 each row
                    % left, right, top, bottom

stPlot.xTick        = 1:nSemitones;
stPlot.yTick        = 0:20:100;
stPlot.yLim         = [0 120];
stPlot.xLim         = [0.5 4.5];
 
% xTicks: frequencies to be displayed along the x-axis
[tones, xTicks]     = get_tones_PR(regs);   stPlot.xTickLabel       = xTicks(:,2:end);

slev        = get_significance_level(nSubjects*nRepetitions, nStrategies); stPlot.slev = slev;

stPlot      = Ensure_field(stPlot,   'SeriesLabel', {'F0mod', 'ACE'});
TitleHead   = ['Pitch Ranking, Ref. '];     stPlot.TitleHead = TitleHead;
stPlot      = Ensure_field(stPlot,'YLabel','% Correct');
XLabel      = 'Comparison tones [Hz]';      stPlot.XLabel = XLabel;
stPlot.ReverseData = 0; % from worst to best score

if length(regs) == 1
    TitleSuffix = {[regs{1} ' (' int2str(xTicks(1,1)) ' [Hz])']};
else
    for i = 1:length(regs)
        TitleSuffix{1,i} = [regs{i} ' (' int2str(xTicks(i,1)) ' [Hz])'];
    end
end

stPlot.TitleSuffix = TitleSuffix;

% f0mData(4,17:20)=NaN;
% aceData(4,17:20)=NaN;

PlotMeans(stPlot, f0mData, aceData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function Scores_correct = quick_PR(filename1)
% function quick_PR(filename1)
%
% Adapted from: thesis_chapter3PRScores
%
% Dependencies:
%   readPRData.m
%
% Examples:
%       close all, clc, figPRScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'CLARINET'})
%
% Programmed by Alejandro, adapted from Matthias' script thesis_chapter3PRScores
%
% Things to Implement:
%       1. Automatic search of Experiments with ACE and F0mod strategies for 
%          all the found subjects (implement this on this script and figMCScores_xPC
%       2. To check significance level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentDirectory = cd;
List_of_files = {   'get_significance_level.m', ...
                    'get_tones_for_pitch_ranking.m', ... %'confidence_interval.m', ... % 
                    'eval_results.m', ...
                    'plot_PR_results.m'};

if isunix == 1
    UserName = 'alejandro';
else
    UserName = 'r0366612';
end

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files,UserName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    instr{1} = '';
end

if nargin == 0
    [filename_part1, filename_part2] = uigetfile('*','PR: Select an APEX experiment result file');
    filename1 = [filename_part2, filename_part1];
end

% regs            = {'G#2', 'C3', 'E3', 'G#3'};
regs                = {'C3', 'G#3'}; % Reference tones to look for...

nRegisters          = length(regs);

disp(['Running script: ' mfilename])

nDifferentTrials  = 4; % To automate (this is constant though)
nStrategies     = 1; % ACE and F0-mod
nSubjects       = 1;
nReferences     = 1;
nSemitones      = 4;

all_res         = zeros( nSubjects*nReferences, nStrategies*nSemitones             );
all_res_spss    = zeros( nSubjects            , nStrategies*nSemitones*nReferences );

j = 1;

for n=1:nReferences

    [all_res(n, 1:nSemitones)] = plot_PR_results(filename1);
    
end

tmp         = all_res( (1:nReferences), 1:nSemitones )';
tmp         = tmp(:);
all_res_spss(1, 1:nReferences*nSemitones) = tmp';
tmp         = all_res((1:nReferences), nSemitones+1:nSemitones)';
tmp         = tmp(:);
all_res_spss(1, nReferences*nSemitones + 1:nReferences*nSemitones) = tmp';
    

% Columns 5x( 1-4 ); 1 means 1 semitone and 4 means 4 semitones
Scores_correct  = all_res_spss(1:nSubjects,     1:nReferences*nSemitones);

p.nRepetitions  = nDifferentTrials;
p.nStrategies   = nStrategies;
p.nSemitones    = nSemitones;
p.nReferences   = nReferences;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end
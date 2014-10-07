function quick_MC(filename1)

% function figMCScores_xPC(dirExpRes, instr)
% 
% Get plots for MCI Experiments. 
% MCI plot (global result) is similar to Figure  3.10, page 64 (but different format)
% Matthias' thesis. 
%
% Examples:
%       figMCScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'UW'})
% Programmed by Alejandro, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'confidenceInterval.m', ...
                    'get_confusion_matrix.m', ...
                    'get_mci_contours.m', ...
                    'plot_MCI_results.m', ...
                    'pool_mci_results.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    [filename_part1, filename_part2] = uigetfile('*',' MCI Select an APEX experiment result file');
    filename1 = [filename_part2, filename_part1];
end 

instr{1} = '';

regs        = {'G#2'};

for i = 1:length(regs)
    note.note   = regs{i}(1:end-1);
    note.octave = str2num(regs{i}(end));
    Freq{i}     = num2str(round(note2freq(note)));
end

% Intervals = {'1 Interval'};

prefix_experiment = {'MC'};

nSubjects = 1;

foundSubjects = 0;

nIntervals      = 1; 
nContours       = 9; % Similar to nSemitones
nSubjects       = 1;

all_res         = zeros(nSubjects*nIntervals, 2*nContours           );
all_res_spss    = zeros(nSubjects           , 2*nContours*nIntervals);
conf            = zeros(nSubjects*nIntervals, 2*nContours           );
nDifferentTrials = 9; % Repetitions per trial

for j=1:nSubjects
        
    for n=1:nIntervals
        
        [   all_res(((j-1)*nIntervals + n),  1:nContours)] = plot_MC_results(filename1);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end

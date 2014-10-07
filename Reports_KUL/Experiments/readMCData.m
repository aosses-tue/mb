function [aceData, f0mData, p] = readMCData(dirExpRes, regs, instr)
% function [aceData, f0mData, p] = readMCData(dirExpRes, regs, instr)
%
% all_res: N x M matrix 
%          N = nSubjects*nIntervals
%          M = nContours*2 (ACE + F0mod)
%
% Each subject occupies nReferences successive rows e.g. row 1 ... 5 for
% subject 1 in case of 5 references (nReferences)
%
% Programmed by Alejandro, adapted from a Matthias' script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    bInstr = 0;
else
    if length(instr) == 0
        bInstr = 0;
    else
        bInstr = 1;
    end
end

if nargin < 2
   regs = {'Gsh2'}; 
end

for i = 1:length(regs)
    note.note   = regs{i}(1:end-1);
    note.octave = str2num(regs{i}(end));
    Freq{i}     = num2str(round(note2freq(note)));
end

Intervals = {'1','2','3','4'};


if nargin == 0
    if isunix == 1
        dirExpRes = '/home/alejandro/Documenten/MM/data/';
    else
        UserName = 'Administrator';
        display(['Windows username being used: ', UserName])
        dirExpRes = ['C:\Documents and Settings\', UserName, '\Desktop\MM\data\'];
    end
end

prefix_experiment = {'MC'};

list = dir([dirExpRes, '*_*']); % considers as subject each folder containing as separator the character '_'

nSubjects = length(list);

foundSubjects = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:length(prefix_experiment) 
    
    for i = 1:nSubjects 
        list_suff = dir([dirExpRes, list(i).name, filesep, '*-', prefix_experiment{j},'*']); % We will look for subfolders containing the experiment Suffix
        disp(['Scanning subject: ' list(i).name])
        if length(list_suff) ~= 0
            foundSubjects = foundSubjects + 1;
            subjects{foundSubjects,1} = list(i).name;
            if length(list_suff) > 1
                for count = 1:length(list_suff)
                    disp(['Type ' num2str(count) ' for selecting: ' list_suff(count).name]);
                end
                idxSelected = input(['More than one experiment found, please enter a number for choosing one of them: '])
                subjects{foundSubjects,2} = list_suff(idxSelected).name;
            else
                subjects{foundSubjects,2} = list_suff.name;
            end
        end
    end
    
end

disp(' ')

for count = 1:length(subjects)
    disp([subjects{count,1} ' in ' subjects{count,2}]);
end
toDelete = input('Type 1 or 0 in a row array for confirming or discarding data: ');

idx2delete = find(toDelete==0);

for count = length(idx2delete):-1:1
    subjects(idx2delete(count),:) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nIntervals      = length(Intervals); 
nContours       = 9; % Similar to nSemitones
nSubjects       = size( subjects, 1 );

disp(['Number of patients found: ' int2str(nSubjects)]);

cd([dirExpRes]);

cur_dir         = pwd;
all_res         = zeros(nSubjects*nIntervals, 2*nContours           );
all_res_spss    = zeros(nSubjects           , 2*nContours*nIntervals);

for j=1:nSubjects
    cd(cur_dir);
    cd([subjects{j, 1} filesep subjects{j, 2}]);
    
    for n=1:nIntervals
        
        InitialsSubject = KeepCapitalLetters(subjects{j,1}); % To obtain the Subject's abbreviation
        if bInstr == 0
            base_filename = [prefix_experiment{1} '_Ref_' Freq{1} '_' 'Int_' Intervals{n} '-'];
        else
            base_filename = [prefix_experiment{1} '_Ref_' Freq{1} '_' 'Int_' Intervals{n} '-' instr '-'];
        end
        
        [   all_res(((j-1)*nIntervals + n),           1:  nContours), ...
            all_res(((j-1)*nIntervals + n), nContours+1:2*nContours), labelContours] = plot_MC_results(base_filename, InitialsSubject, regs{1}, '', 0);
        
    end
    
    tmp = all_res((1:nIntervals) +(j-1)*nIntervals, 1:nContours)';
    tmp = tmp(:);
    all_res_spss(j, 1:nIntervals*nContours) = tmp';                         % 1 row = 1 subject
    tmp = all_res((1:nIntervals) +(j-1)*nIntervals, nContours+1:2*nContours)';
    tmp = tmp(:);
    all_res_spss(j, nIntervals*nContours + 1:2*nIntervals*nContours) = tmp';
    
end

aceData         = all_res_spss(1:nSubjects,              1:  nContours*nIntervals); 
f0mData         = all_res_spss(1:nSubjects, 9*nIntervals+1:2*nContours*nIntervals);
p.subjects      = subjects;
p.labelContours = labelContours;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
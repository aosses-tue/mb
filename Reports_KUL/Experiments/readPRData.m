function [aceData, f0mData, p, stRove] = readPRData(dirExpRes, regs, instr, stSubject)
% function [aceData, f0mData, p, stRove] = readPRData(dirExpRes, regs, instr, stSubject)
%
% Format being read:
%   Test experiment     PR_Ref_131_CLARINET.xml
%   Result file         PR_Ref_131_CLARINET-ACE-AO.apr
%   Where:
%       PR          - Pitch ranking experiment
%       131         - Reference frequency = 131 Hz
%       CLARINET    - Instrument
%       ACE (or F0m)- used Strategy
%       AO          - Subject (Initials)
%
% The script looks for directories inside the dirExpRes folder containing 
% the '_' character in its name (e.g. ..dirExpRes/Alejandro_Osses/)
% The experiments inside each directory should be located in subfolders
% indicating the date and the suffix '-PR', e.g.,
%       ..dirExpRes/Alejandro_Osses/20130415-PR/
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stRove.roveACE = [];
stRove.roveF0m = [];
stRove.roveACEScore = [];
stRove.roveF0mScore = [];

if nargin < 3
    bInstr = 0;
else
    if length(instr) == 0
        bInstr = 0;
    else
        bInstr = 1;
    end
end

if nargin < 4
    stSubject = [];
end

numRegisters = length(regs);
if ~isfield(stSubject,'truncate_data')
    disp([mfilename '.m: There are ' num2str(numRegisters) ' to be loaded. Do you want to truncate any freq. register?'])
    bTruncate = input(['Default = zeros(1,' num2str(numRegisters) '); % define an array with ones for the registers to truncate: '])
    if length(bTruncate) ~= numRegisters
        bTruncate = zeros(1,numRegisters); % normal hearing listeners taken into account
    end
else
    bTruncate = repmat(stSubject.truncate_data,1,numRegisters);
end

if nargin      == 0
    dirExpRes   = '/home/alejandro/Documenten/MM/data/';
end

for i = 1:length(regs)
    note.note   = regs{i}(1:end-1);
    note.octave = str2num(regs{i}(end));
    Freq{i}     = num2str(round(note2freq(note)));
end

prefix_experiment = {'PR'};

disp(' ')
disp('Type ''nh'' for Normal Hearing listeners or ''ci'' for CI recipients.')
% if ~isfield(stSubject,'firstname')
    typeSubject = input('optionally you can type the first name of the subject, e.g. ''ci-Romain'': ');
% else
%     typeSubject = stSubject.firstname;
% end

list            = dir([dirExpRes, typeSubject '*_*']); % considers as subject each folder containing as separator the character '_'
nSubjects       = length(list);
foundSubjects   = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:length(prefix_experiment)
    for i = 1:nSubjects % Now the list is filtered
        list_suff = dir([dirExpRes, list(i).name, filesep, '*-', prefix_experiment{j},'*']); % We will look for subfolders containing the experiment Suffix
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

for count = 1:size(subjects,1)
    disp([subjects{count,1} ' in ' subjects{count,2}]);
end
toDelete = input('Type 1 or 0 in a row array for confirming or discarding data: ');

idx2delete = find(toDelete==0);

for count = length(idx2delete):-1:1
    subjects(idx2delete(count),:) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDifferentTrials  = 4; % To automate (this is constant though)
nStrategies     = 2; % ACE and F0-mod
nSubjects       = size( subjects, 1 );
tones           = get_tones_PR(regs);

disp(['Number of patients found: ' int2str(nSubjects)]);

[nReferences nSemitones] = size(tones);
nSemitones      = nSemitones - 1;

for i=1:nReferences
   tones{i, 1} = strrep(tones{i, 1}, '_', ''); 
end

cd([dirExpRes]);
cur_dir = pwd;

all_res         = zeros( nSubjects*nReferences, nStrategies*nSemitones             );
all_res_spss    = zeros( nSubjects            , nStrategies*nSemitones*nReferences );

for j           = 1:nSubjects
    cd(cur_dir);
    cd([subjects{j, 1} filesep subjects{j, 2}]); % to the subject's folder
    
    tmp_roveACE = [];
    tmp_roveF0m = [];
    tmp_roveACEScore = [];
    tmp_roveF0mScore = [];
    
    for n=1:nReferences
        
        InitialsSubject = KeepCapitalLetters(subjects{j,1}); % To obtain the Subject's abbreviation
        if bInstr == 0
            base_filename = [prefix_experiment{1} '_Ref_' Freq{n} '_'];
        else
            base_filename = [prefix_experiment{1} '_Ref_' Freq{n} '_' instr '-'];
        end
        
        [tmp_ace, tmp_f0m, tmp_Rove] = plot_PR_results(base_filename, InitialsSubject, tones{n, 1}, '', 0, bTruncate(n));
        
        try
            tmp_roveACE = [tmp_roveACE tmp_Rove.roveACE];
            tmp_roveF0m = [tmp_roveF0m tmp_Rove.roveF0m];

            tmp_roveACEScore = [tmp_roveACEScore tmp_Rove.roveACEScore];
            tmp_roveF0mScore = [tmp_roveF0mScore tmp_Rove.roveF0mScore];
        end
        all_res(((j-1)*nReferences + n),            1:  nSemitones) = tmp_ace;
        all_res(((j-1)*nReferences + n), nSemitones+1:2*nSemitones) = tmp_f0m;
    end
    
    try
        stRove.roveACE = [stRove.roveACE; tmp_roveACE];
        stRove.roveF0m = [stRove.roveF0m; tmp_roveF0m];
        stRove.roveACEScore = [stRove.roveACEScore; tmp_roveACEScore];
        stRove.roveF0mScore = [stRove.roveF0mScore; tmp_roveF0mScore];
    end
    
   tmp = all_res((1:nReferences) +(j-1)*nReferences, 1:nSemitones)';
   tmp = tmp(:)/100*nDifferentTrials;
   all_res_spss(j, 1:nReferences*nSemitones) = tmp';
   tmp = all_res((1:nReferences) +(j-1)*nReferences, nSemitones+1:2*nSemitones)';
   tmp = tmp(:)/100*nDifferentTrials;
   all_res_spss(j, nReferences*nSemitones + 1:2*nReferences*nSemitones) = tmp';
    
end

% Columns 5x( 1-4 ); 1 means 1 semitone and 4 means 4 semitones
aceData = all_res_spss(1:nSubjects,                        1:  nReferences*nSemitones)/nDifferentTrials*100;
f0mData = all_res_spss(1:nSubjects, nReferences*nSemitones+1:2*nReferences*nSemitones)/nDifferentTrials*100;

p.nRepetitions  = nDifferentTrials;
p.nStrategies   = nStrategies;
p.nSemitones    = nSemitones;
p.nReferences   = nReferences;
p.InitialsSubject = InitialsSubject;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
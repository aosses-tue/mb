function [LB_ACE_value, LB_F0m_value,  subjects, value_ACE, value_F0m] = readLBData(dirExpRes, test_regs, instr)
% function [LB_ACE_value, LB_F0m_value,  subjects] = readLBData(dirExpRes, regs, instr)
%
% Format being read:
%   Test experiment     LB bla bla
%   Result file         LB_Com_104_UW-p10-F0m-AO.apr
% %   Test experiment     PR_Ref_131_CLARINET.xml
%   Where:
%       LB          - Loudness Balance experiment
%       104         - Test frequency = 104 Hz (Reference always 131 Hz)
%       UW          - Instrument
%       F0m (or ACE)- used Strategy
%       AO          - Subject (Initials)
%
% The script looks for directories inside the dirExpRes folder containing 
% the '_' character in its name (e.g. ..dirExpRes/Alejandro_Osses/)
% The experiments inside each directory should be located in subfolders
% indicating the date and the suffix '-PR', e.g.,
%       ..dirExpRes/Alejandro_Osses/20130611-LB/
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bAllowNaN = 1;
value_ACE = [];
value_F0m = [];

for i = 1:length(test_regs)
    Freq{i}     = num2str(test_regs(i));
end
 
prefix_experiment = {'LB'};

disp(' ')
disp('Type ''nh'' for Normal Hearing listeners or ''ci'' for CI recipients (or nothing).')
typeSubject     = input('optionally you can type the first name of the subject, e.g. ''ci-Romain'': ');
list            = dir([dirExpRes, typeSubject '*_*']); % considers as subject each folder containing as separator the character '_'
nSubjects = length(list);

foundSubjects = 0;

for j = 1:length(prefix_experiment)
    for i = 1:nSubjects % Now the list is filtered
        list_suff = dir([dirExpRes, list(i).name, filesep, '*-', prefix_experiment{j},'*']); % We will look for subfolders containing the experiment Suffix
        if length(list_suff) ~= 0
            for k = 1:length(list_suff)
                if list_suff(k).isdir == 1 % only if the experiments are in a folder
                    foundSubjects = foundSubjects + 1;
                    subjects{foundSubjects,1} = list(i).name;
                    subjects{foundSubjects,2} = list_suff(k).name;
                end
            end
        end
    end
end

% nDifferentTrials  = 4; % To automate (this is constant though)
nStrategies     = 2; % ACE and F0-mod
if exist('subjects','var')
    nSubjects   = size( subjects, 1 );
else
    nSubjects   = 0;
end

disp([mfilename '.m: Number of patients found: ' int2str(nSubjects)]);

nTestTones = 4;

for j           = 1:nSubjects
    
    if j == 1
        cd([dirExpRes]);
        cur_dir = pwd;
    end
    
    cd(cur_dir);
    cd([subjects{j, 1} filesep subjects{j, 2}]); % to the subject's folder
    
    for n=1:nTestTones
        
        InitialsSubject = KeepCapitalLetters(subjects{j,1}); % To obtain the Subject's abbreviation

        base_filename1 = [prefix_experiment{1} '_Com_' Freq{n} '_' instr '-m'];
        base_filename2 = [prefix_experiment{1} '_Com_' Freq{n} '_' instr '-p'];
        
        [LB_value(1), ref_value1] = get_balance_levels(base_filename1, ['ACE-' InitialsSubject]);
        [LB_value(2), ref_value2] = get_balance_levels(base_filename2, ['ACE-' InitialsSubject]);
        
        value_ACE = [value_ACE [LB_value(1);LB_value(2)]];
        
        if bAllowNaN
            if isnan( LB_value(2) )
                LB_value(2) = [];
            end

            if isnan( LB_value(1) )
                LB_value(1) = [];
            end
        end
        
        if ref_value1 == ref_value2
            LB_ACE_value(j,n) = mean( LB_value );
            if length(LB_value) == 1
                disp([mfilename '.m - Subject: ' num2str(j) '.- ' InitialsSubject ', at least one LB value (ACE Strategy) is not complete, are you aware of this?'])
            end
        else
            disp('ACE: Error, press any button')
            pause()
        end
        
        [LB_value(1), ref_value1] = get_balance_levels(base_filename1, ['F0m-' InitialsSubject]);
        [LB_value(2), ref_value2] = get_balance_levels(base_filename2, ['F0m-' InitialsSubject]);
        
        value_F0m = [value_F0m [LB_value(1);LB_value(2)]];
        
        if bAllowNaN
            if isnan( LB_value(2) )
                LB_value(2) = [];
            end

            if isnan( LB_value(1) )
                LB_value(1) = [];
            end
        end
        
        if ref_value1 == ref_value2
            if length(LB_value) == 1
                disp([mfilename '.m - Subject: ' num2str(j) '.- ' InitialsSubject ', at least one LB value (F0mod Strategy) is not complete, are you aware of this?'])
            end
            try 
                load('useCorrectionF0m')
                if n == 1
                    disp([mfilename '.m - Subject: ' num2str(j) '.- ' InitialsSubject ', Correction being applied'])
                end
            catch
                useCorrectionF0m = 0;
            end
            LB_F0m_value(j,n) = mean( LB_value );
            LB_F0m_value(j,n) = LB_F0m_value(j,n) + useCorrectionF0m;
        else
            disp('F0m: Error, press any button')
            pause()
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
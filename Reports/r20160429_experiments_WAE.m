function [score_final score_final_pair] = r20160429_experiments_WAE(dir, files)
% function [score_final score_final_pair] = r20160429_experiments_WAE(dir, files)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 27/04/2016
% Last update on: 27/04/2016 
% Last use on   : 04/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dir = 'D:\Documenten-TUe\06-Lectures-TUe\BEP-2015-2016\03-Odin\Experiment\webaudioevaluationtool-45df85a336d4\saves\';
% files = {'save-ABC-1.xml','save-ABC-2.xml'};
if nargin == 0;
    dir = 'D:\Databases\dir04-Psychoacoustics\WAE\saves-my-results\';
end
if nargin < 1
    % files = {'save-AO-silence-order.xml'};
    files = {'save-AO-silence-order.xml','save-AO-silence-1.xml'};
    % files = {'save.xml'};
end

pairs = [1 2 3;2 4 1; 5 1 2; 1 3 4; 3 5 1; 4 1 5; 4 2 3; 2 3 5; 5 4 2; 3 5 4]; % put the pairs here, manually
minM = min(min(pairs));
maxM = max(max(pairs));
Ms = minM:maxM;

TotalTrials = 0;

[pairs2 Npairs2] = Get_pairwise_combinations(minM,maxM);

M2 = zeros(size(pairs2,1),1);
Pair2found = zeros(size(pairs2,1),1);

for i = 1:length(files)
    resultfile = [dir files{i}];
    
    [score trialID pres_order] = parseWAE_ABC(resultfile);
    
    idx = find(score==1); % find the index of the chosen sound
    chosen = pairs(idx);
    
    % 1. Scores per sound
    for k = Ms
        if i == 1
            M(k,1) =  sum(chosen == k);
        else
            M(k,1) = M(k,1) + sum(chosen == k);
        end
    end
    
    TotalTrials = TotalTrials + length(trialID);
    
    % 2. Scores per pair
    for k = 1:size(score,1)
        
        tmp_score = sort(find(score(k,:) == 0)); % the sound not belonging to the pair was chosen
        for j = 1:size(pairs2,1)
            
            tmp = sort([find( pairs(k,:) == pairs2(j,1) ) find( pairs(k,:) == pairs2(j,2) )]);
            if length(tmp) == 2
                Pair2found(j) = Pair2found(j) + 1; 
                try
                    if tmp == tmp_score
                        M2(j) = M2(j)+1;
                    end
                catch
                    warning('There seems to be a trial without answer')
                end
            end
            disp('')
        end
        
    end
    
end
NperSound = (TotalTrials*3)/(maxM-minM+1);
score_final = [100*M/NperSound Ms'];
score_final_pair = [100*M2./Pair2found   pairs2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

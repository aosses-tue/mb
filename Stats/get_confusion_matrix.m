function [conf, pres, resp, percentage, percentage_total] = get_confusion_matrix(stimuli, pres, resp, l34)
% function [conf, pres, resp, percentage, percentage_total] = get_confusion_matrix(stimuli, pres, resp, l34)
%
% Inputs:
% 	stimuli - labels for each possible sequence
% 	pres	- presented sequence
% 	resp  - sequence answered
%
% For MCI Tests:
%   stimuli =  'R'    'RFA'    'RFL'    'FA'    'FAR'    'FAFL'    'FL'    'FLFA'    'FLR'
%
% Outputs:
%	conf 	- number of correct answers
%	pres	- passed through (same as input)
%	resp 	- passed through (same as input)
%	percentage

% Programmed by Matthias M., revised by AOV
% Original path: /home/alejandro/Documenten/MATLAB/matthias/src/Matlab/Statistics/get_confusion_matrix.m

if nargin == 0
    return;
end

if(length(pres)~= length(resp))
    error('Number of presentations not equal to responses');
end
repetitions = length(pres)/length(stimuli);     % 2

[pres_size_a, pres_size_b] = size(pres);
indices     = get_occurencesAO(pres, stimuli);
% indices = get_occurences(pres);

% conf = zeros(length(stimuli));                % 9 x 9 Matrix (original by Matthias)
conf        = zeros(1,length(stimuli));         % 1 x 9 Matrix (by AOV)
percentage  = zeros(1,length(stimuli));         % 1 x 9 Matrix

if(l34)
    pos = 1;
else
    pos = 0;
end

counter = 1;
for i=1:length(stimuli)
    idx = zeros(repetitions, 1);
    for k=1:pres_size_a*pres_size_b
        if(strmatch(pres{k}, stimuli{i}, 'exact'))
            idx(counter, 1) = k;
            counter = counter + 1;
        end
    end
    
    for j=1:length(idx);
        if idx(j) ~= 0 % Then coincidences were found...
            if strmatch(resp{idx(j)}, stimuli{i}, 'exact');
%               loc_idx = strmatch(resp{idx(j)}, stimuli, 'exact');
                conf(1, i) = conf(1, i) + 1;
            end
        end
    end
end

% percentage = diag(conf)/repetitions*100;
percentage          =      conf ./    (indices)*100;
percentage_total    =  sum(conf) / sum(indices)*100;

if max(max(percentage_total) > 100)
    display(['Warning: apparently the trials are not presented the same amount of times...'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = get_occurences(indices)

for i=1:length(indices)
    indices{i} = isempty(indices{i});
end
indices         = cell2mat(indices);
indices         = find(indices == 0);

function repetitions_per_trial = get_occurencesAO(indices, stimuli)

[sizea, sizeb]  = size(indices);
repetitions_per_trial = zeros(1,length(stimuli));
for i = 1:length(stimuli)
    for j = 1:sizea
        for k = 1:sizeb
            if strmatch(indices{j,k}, stimuli{i}, 'exact');
                repetitions_per_trial(i) = repetitions_per_trial(i) + 1;
            end
        end
    end
end
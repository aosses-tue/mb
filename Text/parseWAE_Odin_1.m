function M = parseWAE_Odin_1(resultfile)
% function M = parseWAE_Odin_1(resultfile)
%
% 1. Description:
%
% 2. Stand-alone example:
%     dir = 'D:\Documenten-TUe\06-Lectures-TUe\BEP-2015-2016\03-Odin\Experiment\webaudioevaluationtool-45df85a336d4\saves\';
%     resultfile = [dir 'save-AO-proc.xml'];
%     M = parseWAE_Odin_1(resultfile);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    dir = 'D:\Documenten-TUe\06-Lectures-TUe\BEP-2015-2016\03-Odin\Experiment\webaudioevaluationtool-45df85a336d4\saves\';
    resultfile = [dir 'save-AO-proc.xml'];
end

pairs = [0 1; 0 2; 0 3; 0 4; 2 1; 3 1; 4 1; 2 3; 2 4; 3 4]; % put the pairs here, manually
minM = min(min(pairs));
maxM = max(max(pairs));

[score trialID pres_order] = parseWAE_Odin_2(resultfile);

M = nan(length(minM:maxM),length(minM:maxM));
for i = minM:maxM
    if i == 1
        disp('')
    end
    i_el = i + 1;
    for j = minM:maxM
        j_el = j + 1;
        if i ~= j
            idx = find(pairs(:,1)==i & pairs(:,2)==j); % first col
            if length(idx) ~= 0 % value found
                col_i = 1;
                col_j = 2;
            else
                idx = find(pairs(:,2)==i & pairs(:,1)==j); 
                col_i = 2;
                col_j = 1;
            end
            if col_i == 1
                M(i_el,j_el) = score(idx,col_i); % score related to 'target' value
            end
            if col_i == 2
                M(i_el,j_el) = score(idx,col_i);
            end
        end
       
        disp('')
    end
    disp('')
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
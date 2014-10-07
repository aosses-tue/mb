function convertPR2SPSS(aceData, f0mData)
% function convertPR2SPSS(aceData, f0mData)
%
% Converts PR outputs to SPSS format
%
% Programmed by Alejandro Osses, ExpORL, KULeuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To prepare input to SPSS:
% 3 within-subject

strat_merged = [aceData f0mData];

disp('3 within subject, C = comparison tone; B = base tone, 0 or 1 = strategy. Each row corresponds to a subject')
disp('   ')
disp('Subject C1_B1_0 C2_B1_0 C3_B1_0 C4_B1_0 C1_B2_0 C2_B2_0 C3_B2_0 C4_B2_0 C1_B3_0 C2_B3_0 C3_B3_0 C4_B3_0 C1_B4_0 C2_B4_0 C3_B4_0 C4_B4_0 C1_B5_0 C2_B5_0 C3_B5_0 C4_B5_0 C1_B1_1 C2_B1_1 C3_B1_1 C4_B1_1 C1_B2_1 C2_B2_1 C3_B2_1 C4_B2_1 C1_B3_1 C2_B3_1 C3_B3_1 C4_B3_1 C1_B4_1 C2_B4_1 C3_B4_1 C4_B4_1 C1_B5_1 C2_B5_1 C3_B5_1 C4_B5_1')

Subjects = [{'S14 ','S16 ','S12 ','S13 ','S11 '}];

for i = 1:size(strat_merged,1)
    tmp_txt = [];
    for j = 1:size(strat_merged,2)
        tmp_txt = [tmp_txt num2str(round(strat_merged(i,j)*100)/100) ' '];
    end
    disp([Subjects{i} tmp_txt])
    % 1. Copy this to a txt file
    % 2. Import into SPSS
end

disp('  ')
disp('  ')

% 2 within-subject
aceDataR = transpose(aceData); aceDataR = aceDataR(:); 
aceDataR = transpose(reshape(aceDataR,4,size(aceData,1)*size(aceData,2)/4));
f0mDataR = transpose(f0mData); f0mDataR = f0mDataR(:);
f0mDataR = transpose(reshape(f0mDataR,4,size(f0mData,1)*size(f0mData,2)/4));

column_ref = [];
for j = 1:size(aceData,1)
    column_ref = [column_ref; [1:5]'];
end

strat_merged = [aceDataR f0mDataR column_ref]; % base tone. 1 means 104, 2 means 131, etc.

strat_ordered = [];
for j=1:5
    idx = find(column_ref==j);
    strat_ordered = [strat_ordered; strat_merged(idx,:)]; 
end
strat_merged = strat_ordered;

disp('2 within subject')
disp('   ')
disp('Subject C1_ACE C2_ACE C3_ACE C4_ACE C1_F0m C2_F0m C3_F0m C4_F0m base_tone')
    
for i = 1:size(strat_merged,1)
    tmp_txt = [];
        
    for j = 1:size(strat_merged,2)
        tmp_txt = [tmp_txt num2str(round(strat_merged(i,j))) ' '];
    end
    
    idx = mod(i,length(Subjects));
    if ~idx; idx = length(Subjects); end;
    
    disp([Subjects{idx} tmp_txt])
    % 1. Copy this to a txt file
    % 2. Import into SPSS
    % variables: C1_0, C2_0, C3_0, C4_0; C1_1, C2_1, C3_1, C4_1; ref_tone
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
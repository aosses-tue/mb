function Create_Piano_main_20150523
% function Create_Piano_main_20150523
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 18/05/2016
% Last update on: 18/05/2016 
% Last use on   : 20/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_skip = 1;
dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '06-Exp-TUe-1-similarity\02-to-be-used-in-AB-comparison' delim];
dir_stimuli = 'Stimuli-main-20160523';

pairs = Get_pairwise_combinations(1,7);

% idx(:,1) = transpose(randperm(Npairs));
% idx(:,2) = transpose(randperm(Npairs));
% idx(:,3) = transpose(randperm(Npairs));
% idx(:,4) = transpose(randperm(Npairs));
% idx(:,5) = transpose(randperm(Npairs));

idx = [      5     6     2     1    11; ...
            20    18    21     2    15; ...
            21    20    10     5     9; ...
             1     9    12     9    18; ...
            16     1    18    14    16; ...
             6     2    19    10    10; ...
            11     3     3    21     6; ...
             4    16    14    15     3; ...
            19    10    20     7     8; ...
             2     4     7    12     4; ...
            15    21    11    13     5; ...
             8     5     5     6    13; ...
            12    12    16    19    17; ...
             3    13     9     8    14; ...
            17    11     8    17     2; ...
            14     8     1    18     1; ...
             9    14    13     3     7; ...
            18     7     4    20    19; ...
            13    15    15    11    20; ...
            10    17     6    16    12; ...
             7    19    17     4    21];
 
Subject_Nr = 1;

for i = 1:size(idx,2)
    
    pairsRand = pairs(idx(:,i),:);
    pairsRandA = pairsRand;
    pairsRandB = pairsRand;
        
    % Alternating order...
    pairsRandA(1:2:21,1) = pairsRand(1:2:21,2);
    pairsRandA(1:2:21,2) = pairsRand(1:2:21,1);
    
    pairsRandB(2:2:21,1) = pairsRand(2:2:21,2);
    pairsRandB(2:2:21,2) = pairsRand(2:2:21,1);
    
    bunch1 = pairsRandA( 1:11,:); % overlap in pair 11
    bunch2 = pairsRandA(11:21,:);
    bunch3 = pairsRandB( 1:11,:); % overlap in pair 11
    bunch4 = pairsRandB(11:21,:);
    
    %%%
    if Subject_Nr < 10
        fname_label = sprintf('S0%.0f',Subject_Nr);
    else
        fname_label = sprintf('S%.0f',Subject_Nr);
    end
    if i == 1
        il_gen_experiments(bunch1,do_skip,fname_label,dir_where,dir_stimuli);
    else
        il_gen_experiments(bunch1,1,fname_label,dir_where,dir_stimuli);
    end
    Subject_Nr = Subject_Nr + 1;
    
    %%%
    if Subject_Nr < 10
        fname_label = sprintf('S0%.0f',Subject_Nr);
    else
        fname_label = sprintf('S%.0f',Subject_Nr);
    end
    if i == 1
        il_gen_experiments(bunch2,do_skip,fname_label,dir_where,dir_stimuli);
    else
        il_gen_experiments(bunch2,1,fname_label,dir_where,dir_stimuli);
    end
    Subject_Nr = Subject_Nr + 1;
    
    %%%
    if Subject_Nr < 10
        fname_label = sprintf('S0%.0f',Subject_Nr);
    else
        fname_label = sprintf('S%.0f',Subject_Nr);
    end
    if i == 1
        il_gen_experiments(bunch3,do_skip,fname_label,dir_where,dir_stimuli);
    else
        il_gen_experiments(bunch3,1,fname_label,dir_where,dir_stimuli);
    end
    
    Subject_Nr = Subject_Nr + 1;
    
    %%%
    if Subject_Nr < 10
        fname_label = sprintf('S0%.0f',Subject_Nr);
    else
        fname_label = sprintf('S%.0f',Subject_Nr);
    end
    if i == 1
        il_gen_experiments(bunch4,do_skip,fname_label,dir_where,dir_stimuli);
    else
        il_gen_experiments(bunch4,1,fname_label,dir_where,dir_stimuli);
    end
    
    Subject_Nr = Subject_Nr + 1;
    
    disp('')
    
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
 
function il_gen_experiments(bunch,do_skip,fname_label,dir_where,dir_stimuli)

opts = [];
opts.dir_where = dir_where;
opts.dirstimuli = dir_stimuli;

opts.RENAME = [fname_label '-A'];
opts.SubjectID = fname_label;
Create_Piano_ICRA_multi_20160518(bunch(1:2,:),do_skip,opts);

opts.RENAME = [fname_label '-B'];
opts.SubjectID = fname_label;
Create_Piano_ICRA_multi_20160518(bunch(3:5,:),do_skip,opts);

opts.RENAME = [fname_label '-C'];
opts.SubjectID = fname_label;
Create_Piano_ICRA_multi_20160518(bunch(6:8,:),do_skip,opts);

opts.RENAME = [fname_label '-D'];
opts.SubjectID = fname_label;
Create_Piano_ICRA_multi_20160518(bunch(9:11,:),do_skip,opts);

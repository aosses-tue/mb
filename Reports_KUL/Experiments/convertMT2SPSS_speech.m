function convertMT2SPSS_speech(p)
% function convertMT2SPSS_speech(p)
%
% Converts Speech Tests outputs to SPSS format
%
% Programmed by Alejandro Osses, ExpORL, KULeuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p.column_score = input('Type 1-4 for using word/phoneme/cons or vowel scores: ')
p.column_score = 9;
Material    = 'MT';

aceData     = p.aceData(:,p.column_score);
f0mData     = p.f0mData(:,p.column_score);
bGroupSNR   = 1;
% To prepare input to SPSS:
% 3 within-subject

SNR = [10 5];

for i = 1:length(SNR)
    idx_SNR = find(p.aceData(:,p.column_SNR)==SNR(i));
    if i==1
        strat_merged = p.aceData(idx_SNR,p.column_subject);
    end
    strat_merged = [strat_merged aceData(idx_SNR) f0mData(idx_SNR)];
end

switch p.column_score
    case 1
        %disp('subject word_ACE_Q word_F0mod_Q word_ACE_10 word_F0mod_10');
        NormFactor = 10;
    case 9
        disp('subject word_ACE_10 word_F0mod_10 word_ACE_5 word_F0mod_5');
        NormFactor = 1; % Stored processed as percentage
    case 3
        %disp('subject consonant_ACE_Q consonant_F0mod_Q consonant_ACE_own_Q consonant_ACE_10 phoneme_F0mod_10 phoneme_ACE_own_10');
        NormFactor = 20;
    case 4
        %disp('subject vowel_ACE_Q vowel_F0mod_Q vowel_ACE_own_Q vowel_ACE_10 vowel_F0mod_10 phoneme_ACE_own_10');
        NormFactor = 10;
end

for i = 1:size(strat_merged,1)
    tmp_txt = [];
    for j = 1:size(strat_merged,2)
        if j==1
            tmp_txt = [tmp_txt num2str(round(strat_merged(i,j)*100)/100) ' '];
        else
            % trans = strat_merged(i,j); % without
            trans = 2*asin(sqrt(strat_merged(i,j)/100))/pi;
            tmp_txt = [tmp_txt num2str(round(trans*100)/100) ' '];
        end
    end
    disp(tmp_txt)
    % 1. Copy this to a txt file
    % 2. Import into SPSS
end

disp('  ')
disp('  ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
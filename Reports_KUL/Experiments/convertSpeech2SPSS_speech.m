function convertSpeech2SPSS_speech(p)
% function convertSpeech2SPSS_speech(p)
%
% Converts Speech Tests outputs to SPSS format. txt output must be manually
% pasted into a *.txt file and if necessary, some characters have to be 
% replaced
% 
% Programmed by Alejandro Osses, ExpORL, KULeuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    disp([mfilename '.m: using score in column ' num2str(p.column_score)])
    Material = 'LT';
catch
    p.column_score = input('Type 1-4 for using word/phoneme/cons or vowel scores: ')
    Material = 'LL';
end

aceData     = p.aceData(:,p.column_score);
aceDataOwn  = p.aceDataOwn(:,p.column_score);
f0mData     = p.f0mData(:,p.column_score);
bGroupSNR   = 1;
% To prepare input to SPSS:
% 3 within-subject

strat_merged = p.aceData(:,p.column_subject);

if isfield(p,'column_SNR')
    strat_merged = [strat_merged p.aceData(:,p.column_SNR)];
end
    
strat_merged = [strat_merged aceData f0mData aceDataOwn];

if strcmp(Material,'LT')
    disp('subject SRT_ACE SRT_F0mod SRT_baseline');
    NormFactor = 100;
else
    switch p.column_score
        case 1
            disp('subject SNR phoneme_ACE phoneme_F0mod phoneme_ACE_own');
            NormFactor = 10;
        case 2
            disp('subject SNR phoneme_ACE phoneme_F0mod phoneme_ACE_own');
            NormFactor = 30;
        case 3
            disp('subject SNR phoneme_ACE phoneme_F0mod phoneme_ACE_own');
            NormFactor = 20;
        case 4
            disp('subject SNR phoneme_ACE phoneme_F0mod phoneme_ACE_own');
            NormFactor = 10;
    end
end
            

for i = 1:size(strat_merged,1)
    tmp_txt = [];
    for j = 1:size(strat_merged,2)
        if j==1 || (isfield(p,'column_SNR') && j == 2)
            tmp_txt = [tmp_txt num2str(round(strat_merged(i,j)*100)/100) ' '];
        else
            tmp_txt = [tmp_txt num2str(round(strat_merged(i,j)*100*100/NormFactor)/100) ' '];
        end
    end
    disp(tmp_txt)
    % 1. Copy this to a txt file
    % 2. Replace . by , if necessary
    % 3. Replace NaN by . if necessary
    % 4. Import into SPSS
end

disp('  ')
disp('  ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
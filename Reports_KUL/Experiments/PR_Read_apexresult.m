function Res = PR_Read_apexresult(List_of_files, nTrials)

% function Res = PR_Read_apexresult(List_of_files)
%
%   Res{k}.p.MatrixAnswer, 3 column matrix:
%           Col 1       Col 2           Col 3           Col 4           Col 5  
%           Total num   Num of correct  Num of wrong    ID Num          Label
%           of trials   answers         answers         of the trial
%
% Dependencies:
%   a3getresults (Apex TB)
%   a3parseresults (Apex TB)
%   get_mci_contours
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    return; % Testing if included in MATLAB Path
end

if nargin < 2
    nTrials = 60; % Default number of trials
end

fontsize    = 10;
Res{1}.Mean = zeros(nTrials,1);
TotalTrials = 0;

for k = 1:length(List_of_files)
    MatrixAnswer        = [];
    r                   = a3getresults(List_of_files{k});
    s                   = a3parseresults(r);
    Number_of_trials    = length(s);
    Correct_answers     = 0;
    Wrong_answers       = 0;
    MatrixAnswer        = zeros(nTrials, 5); % 'Amount_of_trials' different possibilities

    for i = 1:Number_of_trials
        if length( s(1,i).trial ) == 8  % then Trial label is, e.g., trial10
            s(1,i).numTrial                 = str2num( s(1,i).trial(7:8) );
        else % Length = 7               % then Trial label is, e.g., trial3
            s(1,i).numTrial                 = str2num( s(1,i).trial(7) );
        end
        numTrialMod = mod( s(1,i).numTrial, nTrials);
        if numTrialMod == 0
            numTrialMod                     = nTrials;
        end

        switch s(1,i).corrector
            case 'true'
                Correct_answers             = Correct_answers + 1;
                MatrixAnswer(numTrialMod,2) = MatrixAnswer(numTrialMod,2) + 1;
            case 'false'
                Wrong_answers = Wrong_answers + 1;
                MatrixAnswer(numTrialMod,3) = MatrixAnswer(numTrialMod,3) + 1;
        end
        MatrixAnswer(numTrialMod,1) = MatrixAnswer(numTrialMod,1) + 1;
        
        if MatrixAnswer(numTrialMod, 4)     == 0
            MatrixAnswer(numTrialMod,4)     =  s(1,i).numTrial;
            MatrixLabels{numTrialMod}       =  s(1,i).trial;
        end
    end

    Perc                    = Correct_answers/Number_of_trials*100;

    Res{k+1}.p.MatrixAnswer = MatrixAnswer;
    Res{k+1}.p.MatrixLabels = transpose(MatrixLabels);
    Res{1}.Mean             = Res{1}.Mean + Res{k+1}.p.MatrixAnswer(:,2); % Mean of everything
    TotalTrials             = TotalTrials + Number_of_trials;
end
count = length( find( Res{k+1}.p.MatrixAnswer(:,4) ~= 0 ) );
if count ~= nTrials
    nTrials = count;
    display('Some empty rows discarded...')
end

Res{1}.Mean                 = Res{1}.Mean / (TotalTrials/nTrials) * 100;
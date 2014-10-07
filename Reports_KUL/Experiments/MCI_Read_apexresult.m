function Res = MCI_Read_apexresult(List_of_files)

% function Res = Read_apexresult_MCI(List_of_files)
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

nContours = 9;
contours = get_mci_contours(100);
fontsize = 10;
Res{1}.Mean = zeros(9,1);
TotalTrials = 0;

for k = 1:length(List_of_files)
    MatrixAnswer        = [];
    r = a3getresults(List_of_files{k});
    s = a3parseresults(r);
    Number_of_trials    = length(s);
    Correct_answers     = 0;
    Wrong_answers       = 0;
    MatrixAnswer        = zeros(nContours, 3); % 'Amount_of_trials' different possibilities

    for i = 1:Number_of_trials
        if length( s(1,i).trial ) == 8
            s(1,i).numTrial = str2num( s(1,i).trial(7:8) );
        else % Length = 7
            s(1,i).numTrial = str2num( s(1,i).trial(7) );
        end
        numTrialMod = mod( s(1,i).numTrial, 9);
        if numTrialMod == 0
            numTrialMod = nContours;
        end

        switch s(1,i).corrector
            case 'true'
                Correct_answers = Correct_answers + 1;
                MatrixAnswer(numTrialMod,2) = MatrixAnswer(numTrialMod,2) + 1;
            case 'false'
                Wrong_answers = Wrong_answers + 1;
                MatrixAnswer(numTrialMod,3) = MatrixAnswer(numTrialMod,3) + 1;
        end
        MatrixAnswer(numTrialMod,1) = MatrixAnswer(numTrialMod,1) + 1;
    end

    Perc = Correct_answers/Number_of_trials*100;

    Res{k}.p.MatrixAnswer = MatrixAnswer;
    Res{1}.Mean = Res{1}.Mean + Res{k}.p.MatrixAnswer(:,2);
    TotalTrials = TotalTrials + Number_of_trials;
end
Res{1}.Mean = Res{1}.Mean / (TotalTrials/nContours) * 100;


function [res_values res_fields info] = parseapex_VlMatrix(resultfile)
% function [res_values res_fields info] = parseapex_VlMatrix(resultfile)
%
% parses apex apr files (VlMatrix) containing the following information:
%       experiment, startdate, duration, trials and answers
%
%   resultfile    - input file, string
%   info.startdate
%   info.duration - total test time in seconds
%
% % Example: 
% 
% directory = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Romain_Peeters/20131107-MT/';
% cd(directory);
% resultfile = [directory 'simpellijst_6A-SNRp10-F0m-RP.apr'];
% resultfile = 'simpellijst_6A-SNRp10-F0m-RP.apr';
% [res_values res_fields] = parseapex_VlMatrix(resultfile);
%
% Modified by Alejandro Osses, programmed originally by MH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = [];

files = {[resultfile]};

for k = 1:numel(files)
    filename = files{k};
    domnode = xmlread(filename);
    domnode = domnode.getLastChild;
    
    info.experimentfile     = domnode.getAttribute('experiment_file').toCharArray.';
    general                 = domnode.getElementsByTagName('general').item(0);
    interactive             = general.getElementsByTagName('interactive').item(0);
    info.startdate          = general.getElementsByTagName('startdate').item(0).getTextContent.toCharArray.';
    info.duration_s         = general.getElementsByTagName('duration').item(0).getTextContent.toCharArray.';
    
    if ~isnumeric(interactive)
        try
            item1 = interactive.getElementsByTagName('entry').item(0).getAttribute('new_value').toCharArray.';
            tmp_item1 = strsplit(item1,'-');
            info.processing = tmp_item1{1};
            info.subject    = tmp_item1{2};
        end
        try
            item2 = interactive.getElementsByTagName('entry').item(1).getAttribute('new_value').toCharArray.';
            info.SNR = item2;
        end
    end
    % Parse trials
    trials = domnode.getElementsByTagName('trial');
    nTrials = trials.getLength;
    
    i=1; % row
    info.result_words = zeros(nTrials,5);
    
    for l = 0:nTrials - 1
        
        trial = trials.item(l);
        
        trial_id        = char ( trial.getAttribute('id') );
        answer1          = trial.getElementsByTagName('answer'        ).item(0).getTextContent.toCharArray.';
        answer2          = trial.getElementsByTagName('answer'        ).item(1).getTextContent.toCharArray.';
        answer3          = trial.getElementsByTagName('answer'        ).item(2).getTextContent.toCharArray.';
        answer4          = trial.getElementsByTagName('answer'        ).item(3).getTextContent.toCharArray.';
        answer5          = trial.getElementsByTagName('answer'        ).item(4).getTextContent.toCharArray.';
        
        info.result_words(i,:) = [str2num(answer1),str2num(answer2),str2num(answer3),str2num(answer4),str2num(answer5)];
        info.result_sentences(i) = prod(info.result_words(i,:));
        res_values(i,:) = [{trial_id},{answer1},{answer2},{answer3},{answer4},{answer5}];%, {standard}, {stimulus}, {answer}, {correct_answer}, {result}, responsetime];
        i = i+1;
    end
    
        info.score_words = mean(info.result_words(:))*100;
        info.score_sent  = mean(info.result_sentences(:))*100;
        first_half       = info.result_words(1:nTrials/2,:);
        second_half      = info.result_words(nTrials/2+1:nTrials,:);
        info.reliability = [mean(first_half(:))*100 mean(second_half(:)*100)];
        
    res_fields = [{'id'}, {'ans_naam'}, {'ans_ww'}, {'ans_aantal'}, {'ans_kleur'},{'ans_voorwerp'}];
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp([info.experimentfile ': '])
    disp(['Sentence    | word scores:'])
    disp( 'total sent. | words - first half - second half [%] | diff')
    disp(['         ' num2str(info.score_sent) ' | ' Num2Str([info.score_words info.reliability]) '                |' Num2Str(diff(info.reliability))])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

tmp_res_file = strsplit(resultfile,delim);
tmp_res_file = tmp_res_file{end};
disp(['Result file: ' tmp_res_file])
disp(['SNR = ' num2str(info.SNR) ' dB'])
disp(res_fields)

disp([Num2Str(sum(first_half)) ' | ' num2str(sum(sum(first_half)))])
disp([Num2Str(sum(second_half)) ' | ' num2str(sum(sum(second_half)))])
disp('----------------------------')
disp([Num2Str(sum(info.result_words))  ' | ' num2str(sum(sum(info.result_words)))])

Sent_first_half       = sum(info.result_sentences(1:nTrials/2));
Sent_second_half      = sum(info.result_sentences(nTrials/2+1:nTrials));
disp('Sentences first and second half: ')
disp(Num2Str([Sent_first_half Sent_second_half])) 
disp(['Total sentences: ' num2str(sum(info.result_sentences))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
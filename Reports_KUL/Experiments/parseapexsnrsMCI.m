function [res_values res_fields] = parseapexsnrsMCI(resultfile)

% function [res_values res_fields] = parseapexsnrsMCI(resultfile)
%
% Parses apex apr files and writes a cell array with the fields specified during the script.
% These fields are returned in a cell array (res_fields): id, stimulus, 
% answer, correct_answer, answer, result, responsetime. The corresponding 
% values are returned in the cell array resvalues, having as many rows as 
% trials in the APEX experiment.
%
% Tested with Pitch Ranking Experiment:
%       '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/MCI_AC_104_groep_1-results.apr'
%       '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/MCI_AC_104_groep_2-results.apr'
%       '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/MCI_AC_104_groep_3-results.apr'
%
% Modified by AOV from a Michael's script

files = {[resultfile]};

for k = 1:numel(files)
    filename = files{k};
    domnode = xmlread(filename);
    domnode = domnode.getLastChild;
    
    expfile = domnode.getAttribute('experiment_file').toCharArray.';
    general = domnode.getElementsByTagName('general').item(0);
    interactive = general.getElementsByTagName('interactive').item(0);
    startdate = general.getElementsByTagName('startdate').item(0).getTextContent.toCharArray.';
    duration = general.getElementsByTagName('duration').item(0).getTextContent.toCharArray.';
    if ~isnumeric(interactive)
        age = interactive.getElementsByTagName('entry').item(0).getAttribute('new_value').toCharArray.';
        gender = interactive.getElementsByTagName('entry').item(1).getAttribute('new_value').toCharArray.';
    else
        age = '';
        gender = '';
    end
%     results{k + 1} = sprintf('"%s","%s","%s",%s,"%s",%s', filename, expfile, startdate, age, gender, duration);
    
    % Parse trials
    trials = domnode.getElementsByTagName('trial');
    
    i=1; % row
        
    for l = 1:trials.getLength - 1
        
        startdate = general.getElementsByTagName('startdate').item(0).getTextContent.toCharArray.';
        
        trial = trials.item(l);
        
        trial_id        = char ( trial.getAttribute('id') );
        
        stimulus        = trial.getElementsByTagName('stimulus'      ).item(0).getTextContent.toCharArray.'; 
        answer          = trial.getElementsByTagName('answer'        ).item(0).getTextContent.toCharArray.';
        correct_answer  = trial.getElementsByTagName('correctanswer' ).item(0).getTextContent.toCharArray.';
        result          = trial.getElementsByTagName('result'        ).item(0).getTextContent.toCharArray.';
        responsetime    = trial.getElementsByTagName('responsetime'  ).item(0).getTextContent.toCharArray.';
        
        res_values(i,:) = [{trial_id}, {stimulus}, {answer}, {correct_answer}, {result}, responsetime];
        i = i+1;
    end
    res_fields = [{'id'}, {'stimulus'}, {'answer'}, {'correct_answer'}, {'result'}, {'responsetime'}];
    
end
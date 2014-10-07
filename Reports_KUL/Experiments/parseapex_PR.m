function [res_values res_fields info] = parseapex_PR(resultfile)
% function [res_values res_fields info] = parseapex_PR(resultfile)
%
% Parses APEX apr files having some generic fields: experiment, startdate,
% duration, trials, etc.
% An example is given using a Pitch Ranking experiment, but this script
% can easily be adapted for reading any XML-file
%
%   resultfile    - input file, string
%   info.startdate
%   info.duration - total test time in seconds
%
% % Example: 
% 
% resultfile = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Romain_Peeters/20131107-PR/PR_Ref_165_UW_F0m-LB-F0m-RP.apr'
% [res_values res_fields] = parseapex_PR(resultfile);
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
            age = interactive.getElementsByTagName('entry').item(0).getAttribute('new_value').toCharArray.';
        end
        try
            gender = interactive.getElementsByTagName('entry').item(1).getAttribute('new_value').toCharArray.';
        end
    else
        age = '';
        gender = '';
    end
    % results{k + 1} = sprintf('"%s","%s","%s",%s,"%s",%s', filename, expfile, startdate, age, gender, duration);
    
    % Parse trials
    trials = domnode.getElementsByTagName('trial');
    
    i=1; % row
        
    for l = 1:trials.getLength - 1
        
        startdate = general.getElementsByTagName('startdate').item(0).getTextContent.toCharArray.';
        
        trial = trials.item(l);
        
        trial_id        = char ( trial.getAttribute('id') );
        standard        = trial.getElementsByTagName('standard'      ).item(0).getTextContent.toCharArray.'; 
        stimulus        = trial.getElementsByTagName('stimulus'      ).item(0).getTextContent.toCharArray.'; 
        answer          = trial.getElementsByTagName('answer'        ).item(0).getTextContent.toCharArray.';
        correct_answer  = trial.getElementsByTagName('correct_answer').item(0).getTextContent.toCharArray.';
        result          = trial.getElementsByTagName('result'        ).item(0).getTextContent.toCharArray.';
        responsetime    = trial.getElementsByTagName('responsetime'  ).item(0).getTextContent.toCharArray.';
        
        res_values(i,:) = [{trial_id}, {standard}, {stimulus}, {answer}, {correct_answer}, {result}, responsetime];
        i = i+1;
    end
    res_fields = [{'id'}, {'standard'}, {'stimulus'}, {'answer'}, {'correct_answer'}, {'result'}, {'responsetime'}];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

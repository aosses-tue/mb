function [histos, stimuli, buttons, corrector, nEachTrial, nTotalTrial] = eval_results_MC(filenames, disgard, strategies)
% function [histos, stimuli, buttons, corrector, nEachTrial, nTotalTrial] = eval_results_MC(filenames, disgard, strategies)
%
% Evaluates Apex3 results on Musical Instrument Pitch Ranking experiment. 
% The result-filename is put together as follows: 
% MIPR_F0{which_F0, e.g. C3}_results_{strategy, e.g. ACE}-{which_experiment, e.g. 1}.xml
%
% Inputs:
%   path:      Path to result files
%   filenames: Cell-array containing filenames of result-files (order: F0mod, ACE, ACE512, F0ace)
%
% Outputs:
%   nEachTrial: Number of presentations per trial (total trials / 4)
%
% Programmed by Matthias, revised by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
	return;
end

nExperiments = size(filenames, 2); % Number of result files = number of strategies

stimuli     = cell(nExperiments, 1);
buttons     = cell(nExperiments, 1);
correctors  = cell(nExperiments, 1);
histos      = cell(nExperiments, 1);

for i=1:nExperiments
    
    try
        [stimuli{i, 1}, buttons{i, 1}, corrector{i, 1}] = a3cstresults(filenames{i});
    catch
        error('Make sure that the folders under each patient are named properly: ''yyyymmdd-MC''');
    end
    stimuli{i} = stimuli{i, 1}(1+disgard:end);
    for j=1:length(stimuli{i})
        tmp = strread(stimuli{i}{j}, '%s', 'delimiter', '_');
        stimuli{i}{j} = [tmp{1} '_' tmp{end-1}];
    end
    corrector{i, 1} = corrector{i, 1}(1+disgard:end);
    buttons{i, 1} = buttons{i, 1}(1+disgard:end);
    
    % Generates the histogram for each result file:
    histos{i} = process_stimuli_exp(stimuli{i, 1}, corrector{i, 1}, '_', {'stimulus'}); 
    
end

correct_histo   = 0;

nTotalTrial     = length( stimuli{1} );
nEachTrial      =  nTotalTrial / 9;
 
end
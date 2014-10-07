function Processing_APEX_data(PathExperiments)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Processing_APEX_data(PathExperiments)
%
% Reads all the APEX result files (*.apr) in a directory and plots their 
% confusion matrix. Also prints on screen the names of the files being plotted
%
% Dependencies:
%
%   a3getresults    - APEX TB
%   a3parseresults  - APEX TB
%   confusionmatrix - APEX TB
%   plot_pilot_NH_results.m - adapted from Matthias' script plot_thesis_pilot_NH_results.m
%   get_num_total_trials.m
%   get_num_correct_trials.m
%
% Programmed by AOV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'plot_pilot_NH_results.m'};
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    PathExperiments = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/AOV_2013_Tests/20130311/';
end

r = dir( fullfile([PathExperiments '*apr']) );

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('Loaded experiments:')

if length(r) ~= 0
    for i = 1:length(r)
        filename{i} = r(i).name;
        display(['Experiment ' num2str(i) ': ' filename{i}])
    end
else
    display('No *.apr files found...')
    return;
end
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

for j = 1:length(filename)
    short_name  = name2figname(filename{j});

    results     = a3getresults([PathExperiments filename{j}]); % read <processed> section in APEX Experiment
    s           = a3parseresults(results); % get the answers per trial

    stimuli     = arrayfun(@(a) a.correctanswer, s, 'UniformOutput', false); % Read Struct 's' and extract the field 'correct answer'
    responses   = arrayfun(@(a) a.useranswer   , s, 'UniformOutput', false); % Read Struct 's' and extract the field 'useranswer'
    %             arrayfun( FunctionName       , A1, Name          , Value)

    for i = 1:length(stimuli)
        stimuli_mat(i) = str2num( cell2mat( stimuli(i) ) );
    end

    for i = 1:length(responses)
        responses_mat(i) = str2num( cell2mat( responses(i) ) );
    end

    [mat,labels] = confusionmatrix(stimuli, responses); % Stimuli contains the correct answers

    [labels,I]  = sort(labels);
    mat         = mat(I,I);

    % plotconfusion(stimuli_mat-1, responses_mat-1);
    figure
    [fig_a fig_b] = plotconfusion_handle(stimuli_mat-1, responses_mat-1);

    title(['Confusion Matrix for: ' short_name])
    xlabel('Target output - The Reference Tone was: (0) First tone in the trial; (1) Second tone in the trial')
    ylabel('Patient answer')

    legend([fig_b(1) fig_b(2) fig_b(3) fig_b(9)],'Correct','Wrong','Marginal','Final results','Location','NorthEastOutside') 
    % [stimlabels, percent, count] = a3cst2psy(filename, 'stimulus%d')
end

resPR = plot_pilot_NH_results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end
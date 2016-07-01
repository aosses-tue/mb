function [perc tmp_responses proc_names filenames] = quick_constant_multi(directories2check, hoofd_folder)
% function [perc tmp_staircase proc_names filenames] = quick_constant_multi(directories2check, hoofd_folder)
%
% 1. Description:
%       Since 09/03/2015 this script is compatible to extract results from
%       multiprocedure experiments.
%       Valid for 1 directory at a time
% 
% 2. Stand-alone example:
% % 2.1 Example 2016 (Windows computer), using one input file:
%       file = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-4-Multiprocedure\piano-P1-P7-A4-multi-test.apr';
%       quick_constant_multi(file);      
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: r20150313_update.m
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2012-2013
% Created in    : 2012-2013
% Last update on: 30/03/2016 
% Last use on   : 30/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    hoofd_folder = pwd;
end

for j = 1
    
    if isdir(directories2check)
        % files = dir([hoofd_folder directories2check{j} '*.apr*']);
        files = dir([hoofd_folder directories2check '*.apr*']);
    else
        files{1} = directories2check;
    end
    
    for i = 1:length(files)
        
        if isdir(directories2check)
            % file = [hoofd_folder directories2check{j} files(i).name];
            file = [hoofd_folder directories2check files(i).name];
        else
            file = files{1};
        end
        % error('Programming in progress')
        [tmp_responses stimuli procID] = a3constantresults(file);
        Nprocedures = length(fieldnames(tmp_responses));
            
        proc_names = fieldnames(tmp_responses);
        for k = 1:Nprocedures
            values = tmp_responses.(proc_names{k});
            perc(k) = 100*sum(values)/length(values);
            
            fprintf('Procedure %s: %.2f %% (tot. trials = %.0f)\n', proc_names{k}, perc(k), length(values));
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
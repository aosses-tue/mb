function [perc tmp_staircase filenames] = quick_constant_multi(directories2check, opts, dest_folder, hoofd_folder, bSave)
% function [perc tmp_staircase filenames] = quick_constant_multi(directories2check, opts, dest_folder, hoofd_folder, bSave)
%
% 1. Description:
%       Since 09/03/2015 this script is compatible to extract results from
%       multiprocedure experiments.
%       Valid for 1 directory at a time
% 
% 2. Stand-alone example:
% % 2.1. Example 2015:
%       dir2check = {'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\'};
%       quick_staircases(dir2check);
% 
% % 2.2.  Example 2012-2013 (Ubuntu computer):
%       directories2check = {   'ci-Jean-Baptiste_Daumerie/20131016-LT/', ...
%                               'ci-Jean-Baptiste_Daumerie/20131022-LT/', ...
%                               'nh-Anneke_Lenssen/20130806-LT/'};
%       hoofd_folder = '~/Documenten/Meas/Meas/Experiments/Results_XML/';
%       dest_folder  = [hoofd_folder 'ci_pooled/Figures/'];
%       quick_staircases(directories2check, dest_folder, hoofd_folder);
%
% % 2.3 Example 2016 (Windows computer), using one input file:
%       file = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-4-Multiprocedure\piano-P1-P7-A4-multi-test.apr';
%       opts.bPlot = 1; % if you want to plot the figure
%       opts.bSave = 1; % if you want to save the figure
%       quick_staircases(file,opts);      
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

if nargin < 5
    bSave = 0;
end

if nargin < 4
    hoofd_folder = '';
end

if nargin < 3
    dest_folder = Get_TUe_paths('outputs');
end

if nargin < 2
    opts = [];
end

opts = Ensure_field(opts,'mode','mean'); % 'median'
if nargout == 0
    opts = Ensure_field(opts,'bPlot',1);
else
    opts = Ensure_field(opts,'bPlot',0);
end
opts = Ensure_field(opts,'N4mean',6);
opts = Ensure_field(opts,'filter','');

mode    = opts.mode;
N4mean  = opts.N4mean;
bPlot   = opts.bPlot;

N        = 10;

for j = 1
    
    if isdir(directories2check)
        % files = dir([hoofd_folder directories2check{j} '*.apr*']);
        files = dir([hoofd_folder directories2check opts.filter '*.apr*']);
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
        [tmp_staircase stimuli procID] = a3constantresults(file);
        Nprocedures = length(fieldnames(tmp_staircase));
        for k = 1:Nprocedures
            N = max(N, length(tmp_staircase.(procID{k})));
        end
            
        proc_names = fieldnames(tmp_staircase);
        for k = 1:Nprocedures
            values = tmp_staircase.(proc_names{k});
            perc(k) = 100*sum(values)/length(values);
            
            fprintf('Procedure %s: %.2f %% (tot. trials = %.0f)\n', proc_names{k}, perc(k), length(values));
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
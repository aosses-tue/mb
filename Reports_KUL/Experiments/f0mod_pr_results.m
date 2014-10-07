function alldata = f0mod_pr_results(infodir)
% function alldata = f0mod_pr_results(infodir)
%
% Run the following:
%   addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Utility/'); % getpaths.m is here
%   addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Apex_MT/'); % APEX TB
%   addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/Apex_MT/tools/');
%
% Comments added by AO, see following lines:
%       line 34, 35
%       Interesting script: a3cst2psy.m
%
% Programmed by Tom Francart, commented and edited by Alejandro Osses
% Last updated: 12/06/2014 (for Windows-TU/e compatibility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPlot = 1;
bWrite = 0; % set to 1 in case the txt-file does not exist yet in your PC

if nargin == 0
    infodir = [];
end

if isunix
    
    info = getpaths;
    exproot = [info.svn.Meas 'Experiments/Results_XML/'];
    % infodir.export_folder = '~/repo/alejandro/documenten/paper/data/';
    infodir = Ensure_field(infodir, 'export_folder', '~/Documenten/LaTeX_Docs/paper/data_new/'); % local folder Alejandro
    % infodir.export_folder = '~/repo/alejandro/documenten/paper/figures/';
    infodir = Ensure_field(infodir, 'figures_folder', '~/Documenten/LaTeX_Docs/paper/figures/'); % local folder Alejandro

else
    exproot = Get_TUe_paths('ex_APEX_results');
    outputdir = Get_TUe_paths('outputs');
    infodir = Ensure_field(infodir, 'export_folder' , outputdir); 
    infodir = Ensure_field(infodir, 'figures_folder', outputdir);
end


expfolders = {
    'ci-Jan_Leys/20131108-PR'
    'ci-Jean-Baptiste_Daumerie/20131029-PR'
    'ci-Maria_Brughmans/20131125-1007-PR'
    'ci-Patrick_Meul/20131106-PR-retest'
    'ci-Romain_Peeters/20131107-PR'
    'ci-Wouter_David/20131127-PR' };
    
%        ci-Julia_Schoolmeesters/  ? % Julia did not participate on PR experiments
%        ci-Wouter_David/            % 2 incomplete experiments

names = {'subject' 'register' 'semitones' 'strategy' 'percent'};
alldata = dataset([], [], [], [], [], 'VarNames', names);


for f=1:length(expfolders)
    files = dir( [exproot expfolders{f} '/*.apr'] );
    
    for ifile=1:length(files)
        filename = [exproot expfolders{f} '/' files(ifile).name];
        
        names = regexp(filename, '_Ref_(?<register>\d+)_UW_(?<strategy>\w{3})-.*-(?<name>\w{2,3}).apr', 'names');
        if (isempty(names))
            disp(['Skipping file ' filename]);
            continue;
        end
        subject = names.name;
        register = str2num(names.register);
        strategy = names.strategy;
        
        [stimlabels, percent, count, totals] = a3cst2psy(filename, 'stimulus%d');
        
        for i=1:length(stimlabels)
           alldata = [alldata; 
            dataset({subject}, register, stimlabels(i), {strategy}, percent(i), ...
            'VarNames', get(alldata, 'VarNames'))];
        end
    end
end

if bWrite
    export(alldata, ...           % this is file name 
                'file'          , [infodir.export_folder 'CI-PR.txt'], ...
                'Delimiter'     , '\t', ...
                'WriteVarNames' , true);
end

if bPlot % added by AO
    f0mod_pr_plot(alldata, infodir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
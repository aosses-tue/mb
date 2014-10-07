function plottingErrors
% function plottingErrors
%
% Run the script paper_IJA_2014 to see how to use this function
%
% Programmed by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = getpaths('',1);

addpath([info.svn.MATLAB 'Benchs_F0mod_NMT' delim])
addpath([info.svn.MATLAB 'Classes' delim])

info.root_location = '/media/Elements/orl-wrk-0089/Documenten/LaTeX_Docs/predoc_plan_wekelijkse_update/wu028_2013_07_04/Variables/LIST-Err-Auto055/';
info.outputdir = '~/Documenten/LaTeX_Docs/paper/figures/';

figErrValidation(info)
close all

end

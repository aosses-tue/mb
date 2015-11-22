function [y ypd] = Get_filenames(directory, exp2filter, extra)
% function [y ypd] = Get_filenames(directory, exp2filter, extra)
%
% 1. Description:
%   Get file names using a filter specified by exp2filter. The y-variable 
%   corresponds to a cell-array
% 
%   directory     - Directory to be checked
%   exp2filter    - expression to filter file names (e.g. 'sb*.wav', i.e. all
%                   the files starting with 'sb' having wav extension)
%
%   y             - cell array
% 
% 2.1 Example Unix-based:
%   directory        = '~/Documenten/fda_eval_VlMatrix/wav/';
%   extra.bExtension = 0; % To delete extension
%   file_orig        = Get_filenames(directory,['*.wav'],extra);
% 
% 2.2 Example cross-platform:
%   directory = uigetdir('Select a directory with wav-files');
%   file_orig = Get_filenames(directory,['*.wav']);
% 
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2014
% Created in     : 2013-2014
% Last update on : 30/07/2014
% Last use on    : 30/07/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    extra.bExtension = 1;
end

if nargin < 2 
    exp2filter = '';
end

ytmp = dir([directory delim '*' exp2filter]);

if length(ytmp) == 0
    disp('No files found');
    y = [];
end

for i = 1:length(ytmp)
    if extra.bExtension 
        y{i}   = ytmp(i).name;
    else
        y{i} = Delete_extension( ytmp(i).name, exp2filter);
    end
    ypd{i} = [directory delim y{i}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
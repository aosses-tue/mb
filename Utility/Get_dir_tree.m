function pathNames = Get_dir_tree(directory)
% function pathNames = Get_dir_tree(directory)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 29/6/2014
% Last update: 29/6/2014 % Update this date manually
% Last used: 29/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    directory = ['/home/alejandro/Documenten/MATLAB/MATLAB_TUe'];
    % directory = '/home/alejandro/Documenten/MATLAB/MATLAB_TUe/Utility/'; 
end

Count = 0;
CountHidden = 0;

[pathNames, dirNames, fileNames] = dirwalk(directory);

if nargout == 0
    for i = 1:length(pathNames)
        disp(pathNames{i})
        for j = 1:length(fileNames{i})
            
            if ~Has_character(fileNames{i}{j},'~') % Hidden files Windows
                disp(['          -' fileNames{i}{j}])
                Count = Count + 1;
            else
                CountHidden = CountHidden + 1;
            end
            
        end
        
        if length(fileNames{i}) == 0
            disp(['          - Empty directory'])
        end
        
    end
end

fprintf('Total of %.0f files',Count)

if CountHidden ~= 0
    fprintf('\nTotal of %.0f hidden files\n',CountHidden)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
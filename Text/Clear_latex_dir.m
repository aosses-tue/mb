function Clear_latex_dir(dir)
% function Clear_latex_dir(dir)
%
% 1. Description:
%       To be used to clean up LaTeX folders containing compiled files (*.bib,
%       *.aux, etc.). PDFs are also moved to a special folder. 
% 
% 2. Stand-alone example:
%       dir2clean = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-09-16-Vienna-talk\';
%       Clear_latex_dir(dir2clean);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/06/2015
% Last update on: 24/06/2015 % Update this date manually
% Last use on   : 26/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    try
        dir = uigetdir(Get_TUe_paths('lx_Text'));
    catch
        dir = uigetdir(pwd);
    end
    if ~strcmp( dir(end),delim )
        dir = [dir delim];
    end
end

bFirstTimePDF = 1;
bFirstTime    = 1;
files = Get_filenames(dir);

dir1name = [dir 'Cleaned' delim];
dir2name = [dir 'PDFs'    delim];

for i = 1:length(files)
    try
        extension = strsplit(files{i},'.');
        if length(extension) > 1 % then it is cell
            extension = extension{end};
        else
            extension = '';
        end
    end
    try
        if length(extension) ~= 0
            switch extension
                case 'pdf'
                    if bFirstTimePDF
                        bDir2 = Mkdir(dir2name);
                        if bDir2 == 0
                            error(sprintf('Directory %s already exists, rename it or remove it and run this MATLAB script again',dir2name));
                        end
                        bFirstTimePDF = 0;
                    end
                    movefile([dir files{i}],[dir2name files{i}]);
                case 'tex'
                    % Nothing to be done
                case ''
                    % Nothing to be done
                otherwise
                    if bFirstTime
                        bDir1 = Mkdir(dir1name);
                        if bDir1 == 0
                            error(sprintf('Directory %s already exists, rename it or remove it and run this MATLAB script again',dir1name));
                        end
                        bFirstTime = 0;
                    end
                    movefile([dir files{i}],[dir1name files{i}]);
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

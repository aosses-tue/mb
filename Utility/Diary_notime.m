function Diary_notime(name, bDoDiary, directory)
% function Diary_notime(name, bDoDiary, directory)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 01/07/2014
% Last update on: 20/09/2014 % Update this date manually
% Last use on   : 20/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    bDoDiary = 1;
end

if bDoDiary == 0 
    disp(['Diary (log) for function/script ' name ' disabled. Change bDoDiary to 1 to enable this feature...'])
    return;
end

if nargin < 3
    directory = [Get_TUe_paths('outputs') name 'delim'];
end

p = Get_date;

filename = [directory name];

diary(filename)
disp([name ' : log started...'])
disp(['Output file: ' filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

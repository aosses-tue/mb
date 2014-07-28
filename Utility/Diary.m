function Diary(name, directory)
% function Diary(name, directory)
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
% Last update on: 03/07/2014 % Update this date manually
% Last used on  : 03/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = Get_date;

if nargin < 2
    directory = Get_TUe_paths('outputs');
end

filename = [directory 'log-' name '-' p.date2print];
disp(['Don''t forget to set diary to off by the end of the function/script ' name])

diary(filename)
disp([name 'log started...'])
disp(['Output file: ' filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
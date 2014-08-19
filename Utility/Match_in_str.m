function bMatch = Match_in_str(str2match,str_ref)
% function bMatch = Match_in_str(str2match,str_ref)
%
% 1. Description:
%       bMatch will be 1 if all the characters in str2match are present in
%       str_ref
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       str2match = 'TUe';
%       str_ref = 'TUe, een universiteit';
%       Match_in_str(str2match,str_ref) % expected result: bMatch = 1
% 
%       str2match = 'Technische Universiteit Eindhoven';
%       str_ref = 'TUe, een universiteit';
%       Match_in_str(str2match,str_ref) % expected result: bMatch = 0
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 19/08/2014
% Last update on: 19/08/2014 % Update this date manually
% Last use on   : 19/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length( strfind(str_ref,str2match) )~= 0
    bMatch = 1;
else
    bMatch = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

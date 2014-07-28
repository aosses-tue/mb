function bHas = Has_character(name,Char1)
% function bHas = Has_character(name,Char1)
%
% 1. Description:
%       Looks for Char1 along the 'name' string returning 1 if Char1 is found
%       and 0 otherwise.
% 
% 2. Additional info:
%   Tested cross-platform:
%       Ubuntu: YES, on 17/07/2014
%
% 3. Stand-alone example:
%       txt = '~casa';
%       bHas = Has_character(txt,'~');
% % bHas will be 1, since the string '~casa' it does contain a '~' char
%   
%       txt = '~casa';
%       bHas = Has_character(txt,'/');
% % bHas will be 0, since the string '~casa' it does not contain a '/' char
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 17/07/2014
% Last update on: 17/07/2014 % Update this date manually
% Last used on  : 17/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_name    = name;
cont        = find(name == Char1);
bHas        = 0;

if cont ~= 0
    bHas = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
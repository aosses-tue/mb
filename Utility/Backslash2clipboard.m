function Backslash2clipboard
% function Backslash2clipboard
%
% 1. Description:
%       Copies backslash ('\') to clipboard. This is useful when using an
%       AZERTY keyboard layout that does not contain the backslash.
% 
% 2. Additional info:
%   Tested cross-platform: 
%       Ubuntu: YES
%
% 3. Stand-alone example:
%       Backslash2clipboard;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 17/07/2014
% Last update on: 17/07/2014 % Update this date manually
% Last used on  : 17/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clipboard('copy','\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
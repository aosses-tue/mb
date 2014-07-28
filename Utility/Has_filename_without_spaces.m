function [y bCompatible] = Has_filename_without_spaces(x)
% function [y bCompatible] = Has_filename_without_spaces(x)
%
% 1. Description:µ
%   Detects whether a string contains or not a space in its file name, this 
%   is not-compatible with some softwares
%
%   x - corresponds to a character (letter, words or sentences)
%   y - compatible file name (without spaces)
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 19/5/2014
% Last update: 19/5/2014 % Update this date manually
% Last used: 19/5/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bCompatible = 1;

for i = 1:length(x)
    if strcmp(x(i),' ')
        bCompatible = 0;
        y(i) =  '-';
    else
        y(i) =  x(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
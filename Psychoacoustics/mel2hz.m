function f = mel2hz(m)
% function f = mel2hz(m)
%
% 1. Description:
%       http://code.soundsoftware.ac.uk/projects/ishara/embedded/dsp/mel2hz.html
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/08/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=700*((10.^(m/2595))-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function m = hz2mel(f)
% function m = hz2mel(f)
%
% 1. Description:
%       http://code.soundsoftware.ac.uk/projects/ishara/embedded/dsp/hz2mel.html
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

m=2595*log10(1+f/700);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

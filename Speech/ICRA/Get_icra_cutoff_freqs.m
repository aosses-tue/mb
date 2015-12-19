function cross = Get_icra_cutoff_freqs(type)
% function cross = Get_icra_cutoff_freqs(type)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 17/12/2015
% Last update on: 17/12/2015 
% Last use on   : 17/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    type = 0;
end

switch type
    case 0 % As in ICRA for speech
        cross   = [800 2400]; % same as used in il_getfilters
    case 1
        cross   = [400 800 2400 4800];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

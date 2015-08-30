function param = Get_envelope_descriptors(insig,type)
% function param = Get_envelope_descriptors(insig,type)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       See also LNN_noise_NW.m (where these descriptors have been utilised)
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 28/08/2015
% Last update on: 28/08/2015 
% Last use on   : 28/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    type = 'W';
end
env=abs(hilbert(insig));

switch type
    case 'W'
        param = mean(insig.^4)/( mean(insig.^2).^2 );
    case 'V'
        param = 20*log10(std(env)/mean(env));
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

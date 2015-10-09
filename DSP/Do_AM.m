function [outsig env] = Do_AM(insig,fs,fmod,mdepth)
% function [outsig env] = Do_AM(insig,fs,fmod,mdepth)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/10/2015
% Last update on: 06/10/2015 
% Last use on   : 06/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dBFS = 100;

N = length(insig);
% env     = zeros(N,1); % Initialisation
wstep	= 2*pi/fs;

t = wstep*fmod*(1:N);

env = transpose(1+(mdepth*sin(t)));
outsig = insig.*env;

lvl = rmsdb(insig)+dBFS;
outsig = setdbspl(outsig,lvl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

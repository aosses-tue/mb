function y = dw_rms(x)
% function y = dw_rms(x)
% 
%  1. Description:
%       RMS value of array x
% 
% This file belongs to the roughness algorithm
% contact for the original source code :
% http://home.tm.tue.nl/dhermes/

% Included into psysound by Matt Flax <flatmax @ http://www.flatmax.org> : Matt Flax is flatmax
% March 2007 : For the psySoundPro project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m   = length(x);
y   = 0;

try
    y = sqrt(1/m*sum(x.*x));
catch
    warning('Using code in its original ''not optimised'' form...');
    for k=1:1:m
       y=y+(x(k)^2);
    end
    y = sqrt(y/m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function [y centres] = Probability_density_function(x,N)
% function [y centres] = Probability_density_function(x,N)
%
% 1. Description:
%   Probability density function of 'x'
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%   N = 20;
%   x = wgn(100,1,1);
%   y = Probability_density_function(x,N);
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 21/5/2014
% Last update: 21/5/2014 % Update this date manually
% Last used: 21/5/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    N = 10;
end

[p centres] = hist(x,N);

y = p/sum(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

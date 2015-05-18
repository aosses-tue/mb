function [y centres] = Probability_density_function(x,N)
% function [y centres] = Probability_density_function(x,N)
%
% 1. Description:
%   Probability density function of 'x'. The factor cal ensures that the 
%   probabilioty density function has an area under the curve equal to 1.
%
% 2. Stand-alone example:
%   N = 20;
%   x = wgn(100,1,1);
%   y = Probability_density_function(x,N);
%
% 3. Additional info:
%   Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 21/05/2014
% Last update on: 16/05/2015 % Update this date manually
% Last use on   : 16/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    N = 50;
end

[p centres] = hist(x,N);

dStep   = centres(2)-centres(1);
cal     = 1/(dStep * sum(p));

y = cal * p;

if nargout == 0
    
    figure;
    plot(centres,y); grid on
    title('Probability density function')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

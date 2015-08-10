function [y centres stats] = Probability_density_function(x,N)
% function [y centres stats] = Probability_density_function(x,N)
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
% Last update on: 16/05/2015 
% Last use on   : 10/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    N = 50;
end

if length(N) == 1
    [p centres] = hist(x,N);
else
    centres = N;
    p = hist(x,centres);
end

dStep   = centres(2)-centres(1);
cal     = 1/(dStep * sum(p));

y = cal * p;

if nargout > 2
    idx = min(find(y == max(y))); % in case two maximum values are found, the minimum index (arbitrary) is chosen
    stats.mean  = mean(x);
    stats.std   = std(x);
    
    if kstest( (x-mean(x))/std(x) )
        stats.bNormal = 'no';
    else
        stats.bNormal = 'yes';
    end
end

if nargout == 0
    
    figure;
    plot(centres,y); grid on
    title('Probability density function')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

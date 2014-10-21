function y = Normalise_signal(x)
% function y = Normalise_signal(x)
%
% 1. Description:
%       Normalisation of a signal x. If x has more than 1 column, each 
%       column is interpreted as a different time series
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       x = wgn(1,10,1); % 10 element white noise
%       y = Normalise_signal(x);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/10/2014
% Last update on: 20/10/2014 % Update this date manually
% Last use on   : 20/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = nan(size(x)); % just allocation

if size(x,1) ~= 1
    for i = 1:size(x,2)
        y(:,i) = normalize(x(:,i),'energy');
    end
else
    % in case of 1 row
    y = normalize(x,'energy');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

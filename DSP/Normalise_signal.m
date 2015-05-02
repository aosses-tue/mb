function y = Normalise_signal(x,fs)
% function y = Normalise_signal(x,fs)
%
% 1. Description:
%       Normalisation of a signal x. If x has more than 1 column, each 
%       column is interpreted as a different time series.
%       This normalisation is approximately independent respect to the 
%       sampling frequency.
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       x = wgn(1,10,1); % 10 element white noise
%       fs = 44100;
%       y = Normalise_signal(x,fs);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 20/10/2014
% Last update on: 30/04/2015 % Update this date manually
% Last use on   : 30/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = nan(size(x)); % just allocation

DeltaT = 1/fs;
if size(x,1) ~= 1
    for i = 1:size(x,2)
        CalFactor = sqrt( 1/(DeltaT* sum(x(:,i).^2) ));
        y(:,i) = x(:,i) * CalFactor;
    end
else
    CalFactor = sqrt( 1/(DeltaT* sum(x.^2) ));
    % in case of 1 row
    y = x * CalFactor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

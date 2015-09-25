function [y y4optdet y3] = Normalise_signal(x,fs,method)
% function [y y4optdet y3] = Normalise_signal(x,fs,method)
%
% 1. Description:
%       Normalisation of a signal x. If x has more than 1 column, each 
%       column is interpreted as a different time series.
%       This normalisation is approximately independent respect to the 
%       sampling frequency.
% 
%       This script makes the multiplication sum(y^2) equal to (N/fs), where
%       N/fs is the duration of the digital signal in seconds.
% 
%       If x (at fs) has a duration of 3.65 [s] then sum(y.*y) will be 3.65
%       For the same signal x, optimaldetector(y4optdet,y4optdet) = 1
% 
% 2. Stand-alone example:
%       x = wgn(1,10,1); % 10 element white noise
%       fs = 44100;
%       y = Normalise_signal(x,fs);
% 
%       fs = 44100;
%       x = wgn(1,3.65*fs,1); % 3.65-s length white noise
%       [y y4opt] = Normalise_signal(x,fs);
%       disp(num2str(sum(y.*y))) % this should give 3.65 of total energy
%       disp(num2str(optimaldetector(y4opt,y4opt))) % correlation with itself (autocorrelation) = 1
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also optimaldetector.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 20/10/2014
% Last update on: 30/04/2015 
% Last use on   : 24/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    method = 1;
end

y = nan(size(x)); % just allocation

N = length(x);
if size(x,1) ~= 1
    for i = 1:size(x,2)
        CalFactor = sqrt( N/(fs* sum(x(:,i).^2)) );
        y(:,i) = x(:,i) * CalFactor;
        y4optdet(:,i) = x(:,i) * CalFactor * fs^(1/4);
        switch method
            case 1
                c = sqrt( 1/sum(x(:,i).^2) );
            case 2
                c = sqrt(fs/sum(x(:,i).^2) );
        end
        y3(:,i) = c*x(:,i);
    end
else
    CalFactor= 1/sqrt(fs) * sqrt( N/(sum(x.^2)));
    y        = x * CalFactor;
    y4optdet = x * CalFactor*fs^(1/4);
    switch method
        case 1
            c = sqrt( 1/sum(x.*x) );
        case 2
            c = sqrt(fs/sum(x.*x) );
    end
    y3 = c*x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

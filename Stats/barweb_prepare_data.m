function [mean2plot, std2plot] = barweb_prepare_data(x,y)
% function [mean2plot, std2plot] = barweb_prepare_data(Data)
%
% 1. Description:
%       Calculates mean and std for a set of x- and y-data, assumes that x
%       and y have the same dimensions and also that rows contain different
%       samples/trials for the same condition, calculating M averages and 
%       standard deviations, with M = size(x,1).
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       [m s] = barweb_prepare_data(x,y);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 2/6/2014
% Last update: 2/6/2014 % Update this date manually
% Last used: 2/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = transpose(x);
if nargin == 2
    y = transpose(y);
end

mx = mean(x);
sx = std(x);
if nargin == 2
    my = mean(y);
    sy = std(y);
end

if nargin == 2
    mean2plot   = [mx; my];
    std2plot    = [sx; sy];
else
    mean2plot   = [mx];
    std2plot    = [sx];
end
mean2plot   = transpose(mean2plot);
std2plot    = transpose(std2plot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
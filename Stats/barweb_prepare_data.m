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
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 2/6/2014
% Last update: 2/6/2014 % Update this date manually
% Last used: 2/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = transpose(x);
y = transpose(y);

mx = mean(x);
sx = std(x);
my = mean(y);
sy = std(y);

mean2plot   = [mx; my];
std2plot    = [sx; sy];

mean2plot   = transpose(mean2plot);
std2plot    = transpose(std2plot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
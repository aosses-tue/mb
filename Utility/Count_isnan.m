function [N_is_n N_is_nan] = Count_isnan(x)
% function [N_is_n N_is_nan] = Count_isnan(x)
%
% 1. Description:
%   Counts across a row how many numbers (N_is_n) and how many not a number
%   (N_is_nan) there are...
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 19/6/2014
% Last update: 19/6/2014 % Update this date manually
% Last used: 19/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_is_n      =  transpose(sum( ~isnan(x') ));
N_is_nan    =  transpose(sum(  isnan(x') ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
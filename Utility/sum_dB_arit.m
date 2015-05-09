function Y = sum_dB_arit(X)
% function Y = sum_dB_arit(X)
%
% 1. Description:
%       Y [dB] (related to a pressure of y [Pa]) is obtained by combining 
%       Xn SPL levels (with a pressure of xn [Pa]) in such a way that y is:
%           y = sum(x)
% 
%       see also sum_db.m
% 
% 2. Stand-alone example:
%       X = [60 60];
%       Y = sum_dB_arit(X); % expected value: 66 dB
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/05/2015
% Last update on: 06/05/2015 % Update this date manually
% Last use on   : 06/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x2      = 10.^(X/20); % Anti-log
Y       =  20*log10( sum(x2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

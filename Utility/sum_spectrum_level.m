function [Y ymean] = sum_spectrum_level(X,Delta_f)
% function [Y ymean] = sum_spectrum_level(X,Delta_f)
%
% 1. Description:
%       Y [dB] (related to a pressure of y [Pa]) is obtained by combining 
%       Xn SPL levels (with a pressure of xn [Pa]) in such a way that y is:
%           y = sqrt(sum(x.^2))
% 
%       see also sum_dB_arit.m
% 
% 2. Stand-alone example:
%       X = [60 60];
%       Df = 1; % 1 Hz
%       Y = sum_spectrum_level(X,Df); % expected result: 63 dB
% 
%       X = [60 60 60 60];
%       Df = 0.5; % 0.5 Hz
%       Y = sum_spectrum_level(X,Df); % expected result: 63 dB
%
%       X = 40;
%       Df = 5000; % Hz
%       Y = sum_spectrum_level(X,Df);
% 
% 3. Additional info:
%   Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/06/2014
% Last update on: 06/05/2015 
% Last use on   : 17/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    Delta_f = 1;
end

x2 = 10.^(X/10); % Anti-log squared

Y  =  10*log10( sum(x2).*Delta_f );

if nargout > 1
        
    K = size(X,2);
    for i = 1:K
        idx = find(X(:,i)>floor_dB);
        N = length(idx);
        ymean(i) =  10*log10(1/N * sum(x2(:,i)) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
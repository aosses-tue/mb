function [Y ymean] = sum_db(X,floor_dB)
% function [Y ymean] = sum_db(X,floor_dB)
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
%       Y = sum_db(X); % expected result: 63 dB
% 
% 3. Additional info:
%   Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/06/2014
% Last update on: 06/05/2015 
% Last use on   : 17/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x2 = 10.^(X/10); % Anti-log squared

Y  =  10*log10( sum(x2) );

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
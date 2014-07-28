function [fp Hp] = Get_predominant_values(f,H)
% function [fp Hp] = Get_predominant_values(f,H)
%
% 1. Description:
%       H in dB
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 24/7/2014
% Last update on: 24/7/2014 % Update this date manually
% Last used on  : 24/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b] = size(H);
c     = length(f);

if a == c
    m = b;
    f = f(:); % column vector
else
    m = a;
end

idx = [];
Hp = [];
for i = 1:m
    idx = [idx; find( H(:,i)==max(H(:,i)) )];
    Hp = [Hp; H( idx(i),i)];
end

fp =  f(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function [fp Hp] = Get_predominant_values(f,H,mode)
% function [fp Hp] = Get_predominant_values(f,H,mode)
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
% Created on    : 24/07/2014
% Last update on: 24/07/2014 % Update this date manually
% Last use on   : 30/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    mode = 1; % 'centroid' mode
end

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

switch mode
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0 % old mode, just looking at maximum STFT value
        for i = 1:m
            idx = [idx; find( H(:,i)==max(H(:,i)) )];
            Hp = [Hp; H( idx(i),i)];
        end
        fp =  f(idx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        fp = [];
        for i = 1:m
            %idx = [idx; find( H(:,i)==max(H(:,i)) )];
            fc = Get_spectral_centroid_dB(f, H(:,i),20);
            fp = [fp; fc];
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [y ymean] = sum_db(x,floor_dB)
% function [y ymean] = sum_db(x,floor_dB)
%
% 1. Description:
%   dB sum
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 06/06/2014
% Last update on: 06/06/2014 % Update this date manually
% Last use on   : 10/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x2 = 10.^(x/10); % Anti-log squared

y       =  10*log10( sum(x2) );

if nargout > 1
        
    K = size(x,2);
    for i = 1:K
        idx = find(x(:,i)>floor_dB);
        N = length(idx);
        ymean(i) =  10*log10(1/N * sum(x2(:,i)) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
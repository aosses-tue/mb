function Thres = demo_dau1996b_decision(mue,SPL,criterion)
% function Thres = demo_dau1996b_decision(mue,SPL,criterion)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 18/03/2015
% Last update on: 18/03/2015 % Update this date manually
% Last use on   : 18/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Thres = [];

[N M] = size(mue);

for i = 1:N
    
    Thres(i) = interp1(mue(i,:),SPL,criterion,'linear','extrap');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

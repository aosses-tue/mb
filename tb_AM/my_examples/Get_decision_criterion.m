function out = Get_decision_criterion(ir1,ir2,param,idx)
% function out = Get_decision_criterion(ir1,ir2,param,idx)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 30/04/2015
% Last update on: 30/04/2015 % Update this date manually
% Last use on   : 30/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    idx = 1:size(ir1,2);
end

switch param
    case 'dprime'
        out = mean(ir2(:,idx))-mean(ir1(:,idx));
    case 'cross-correlation'
        out = optimaldetector(ir2(:,idx),ir1(:,idx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

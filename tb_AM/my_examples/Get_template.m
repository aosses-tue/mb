function out = Get_template(ir1,ir2,param,idx)
% function out = Get_template(ir1,ir2,param,idx)
%
% 1. Description:
%       ir - internal representations
% 
% 2. Stand-alone example:
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 23/10/2014 % Update this date manually
% Last use on   : 23/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    idx = 1:size(ir1,2);
end

switch param
    case 'template'
        disp('Make sure ir2 is N+Suprathreshold signal')
        out = Normalise_signal(ir2(:,idx)-ir1(:,idx));
    case 'dprime'
        out = mean(ir2(:,idx))-mean(ir1(:,idx));
    case 'cross-correlation'
        disp('Make sure ir1 is the template')
        out = optimaldetector(ir1(:,idx),ir2(:,idx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

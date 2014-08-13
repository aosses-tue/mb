function y = Get_MIR_param(a,Param)
% function y = Get_MIR_param(a,Param)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       a = miraudio('ragtime');
%       x   = Get_MIR_param(a,'Data');
%       t   = Get_MIR_param(a,'Time');
%       fs  = Get_MIR_param(a,'Sampling');
%       figure;
%       plot(t,x), grid on
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/08/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = get(a,Param);

while iscell(y)
    y = y{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

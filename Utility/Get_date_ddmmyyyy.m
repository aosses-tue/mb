function y = Get_date_ddmmyyyy(b_with_time)
% function y = Get_date_ddmmyyyy(b_with_time)
%
% 1. Description:
%       Returns current 'date' or 'date and time'
% 
% 2. Stand-alone example:
%       Get_date_ddmmyyyy;
% 
% 3. Additional info:
%   Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on : 30/06/2014
% Last update: 30/06/2014 % Update this date manually
% Last used  : 25/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    b_with_time = 0;
end

p = Get_date;

if b_with_time == 0
    y = [p.dd '/' p.mm '/' p.yyyy];
else
    y = [p.dd '/' p.mm '/' p.yyyy ' at ' p.time];
end

if nargout == 0
    disp([p.dd '/' p.mm '/' p.yyyy ' at ' p.time]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
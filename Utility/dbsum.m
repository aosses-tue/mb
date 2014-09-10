function y = dbsum(SPLs)
% function y = dbsum(SPLs)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       SPLs = [60 61];
%       dbsum(SPLs)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/09/2014
% Last update on: 05/09/2014 % Update this date manually
% Last use on   : 05/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = 10*log10(sum(10.^(SPLs/10)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

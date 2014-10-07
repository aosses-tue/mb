function y = dbmean(SPLs)
% function y = dbmean(SPLs)
%
% 1. Description:
%       Energetic mean. Different rows are different times, while different
%       columns refer to different frequency bands
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/10/2014
% Last update on: 03/10/2014 % Update this date manually
% Last use on   : 03/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b] = size(SPLs);

y = 10*log10(1/a*sum(10.^(SPLs/10)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

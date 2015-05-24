function [bNormal y] = Is_normal_distributed(x)
% function [bNormal y] = Is_normal_distributed(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 21/05/2015
% Last update on: 21/05/2015 % Update this date manually
% Last use on   : 21/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bNormal = ~kstest( (x-mean(x))/std(x) );

if bNormal
    y = 'yes';
else
    y = 'no';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

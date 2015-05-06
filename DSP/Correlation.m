function cc = Correlation(x,y,method)
% function cc = Correlation(x,y,method)
%
% 1. Description:
%       If method is 'coefficient', then the correlation coefficient (value
%       between -1 and 1 is returned.
% 
% 2. Stand-alone example:
%       data1 = wgn(1,100,1); % 100 element white noise
%       data2 = wgn(1,100,1); % 100 element white noise
%       CF = Correlation(data1,data2,'coefficient');

% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 01/05/2015
% Last update on: 02/05/2015 % Update this date manually
% Last use on   : 02/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    method = 'coefficient';
end

switch method
    case 'coefficient'
        cfac =	cov(x,y);
        den	=	diag(cfac);
        den	=	sqrt(den*den');
        if den(2,1)>0
          cc =	cfac(2,1)/den(2,1);
        else
          cc =	0;
        end
    otherwise
        error('Correlation method not recognised')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

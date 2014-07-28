function [ydB, num, den] = Get_filter_specs(Hd)
% function [ydB, num, den] = Get_filter_specs(Hd)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       % Run OB_analysis to see an example
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 06/06/2014
% Last update on: 06/06/2014 % Update this date manually
% Last used on  : 17/06/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N       = 4096;
K       = N/2;
h       = [];

[q, num, den, sv, issvnoteq2one] = dispatchsecfilter(Hd);

y = ones(K,1);
for i=1:size(num,1)
    [aux] = freqz(num(i,:),den(i,:),K);
    h = [h aux];
    y = y.*h(:,i);
end

% [y,zf] = df2sosfilter(q,num,den,sv,issvnoteq2one,x,zi);

% ydB = 20*log10(abs(y));
% plot(ydB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
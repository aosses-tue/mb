% Get Gammatone filter coef for DRNL implementation
% Morten L�ve Jepsen
function [b,a]=coefGtDRNL(fc,BW,n,fs);

theta = 2*pi*fc/fs;
phi   = 2*pi*BW/fs;
alpha = -exp(-phi)*cos(theta);

b1 = 2*alpha;
b2 = exp(-2*phi);
a0 = abs( (1+b1*cos(theta)-i*b1*sin(theta)+b2*cos(2*theta)-i*b2*sin(2*theta)) / (1+alpha*cos(theta)-i*alpha*sin(theta))  );
a1 = alpha*a0;

% adapt to matlab filter terminology
b=[a0, a1];
a=[1, b1, b2];

% SE 23.08.2010 15:42 seems to be removed since it was never used with n != 1

% b2=1;
% a2=1;
%
% for j=1:n
% 	b2=conv(b,b2);
% 	a2=conv(a,a2);
% end
%
% b=b2;
% a=a2;

%eof
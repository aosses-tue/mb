function [a0, a1] = LinRegress(x, y)
% LinRegress(x, y)
%     Calculates the linear regression line 'x*a1*x+a0' 
%     through the points defined by 'x' and 'y'.  
%     If division by zero occurs, a1 and a0 become NaN.
[l1, l2] = size(x);
l = l1*l2;
sx  = sum(sum(x));
sy  = sum(sum(y));
sxy = sum(sum(x.*y));
sx2 = sum(sum(x.*x));
% sy2 = sum(sum(y.*y));
d   = l*sx2 - sx.*sx;
if ~(d==0),
    a1 = (l*sxy-sx.*sy)/d;
    a0  = (sy-a1*sx)/l;
else
    a1 = NaN;
    a0 = NaN;
end

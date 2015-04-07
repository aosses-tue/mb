function y = shiftcov(f1, f2, amount)
% Author: V.J.P. Jourdes
% This function calculates the maximum of the cross correlation between two
% signal calculated after shifting one of them of a certain number of
% samples.
%
% y = shifcov(f1, f2, amount);
%
% f1: signal 1
% f2: signal 2
% amount: defines the maximum number of samples that the signal will be
% shifted.

L = length(f1);
ff2 = f2;
x=1;
for i = 1:10:amount, % The signal is shifted 10 by 10 samples for less
    % calculation complexity
    cfac = cov(f1, f2);
    den = diag(cfac);
    den = sqrt(den * den’);
    if den(2,1) > 0
        r(x) = cfac(2,1) / den(2,1);
    else
        r(x) = 0;
    end
    x=x+1;
    f2(1) = [];
    f2(L) = 0;
end

for i = 1:10:amount,
    f1 = [zeros(1,i) f1(1:L-i)];
    cfac = cov(f1, ff2);
    den = diag(cfac);
    den = sqrt(den * den’);
    if den(2,1) > 0
        r(x) = cfac(2,1) / den(2,1);
    else
        r(x) = 0;
    end
    x=x+1;
end
y = max(r);

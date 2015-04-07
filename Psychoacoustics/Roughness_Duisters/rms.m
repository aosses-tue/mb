% out = rms(in);
%
% Calculates the rms value of array 'in'
%
function out = rms(in)

[m, n] = size(in);
out = sqrt(sum(x.^2, 2) / n);
% out = rmsDik(in);
%
% Calculates the rms value of array 'in'
%
function out = rmsDik(in)

[m, n] = size(in);
out = sqrt(sum(in.^2, 2) / n);

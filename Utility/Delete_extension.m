function y = Delete_extension(x, extension)
% function y = Delete_extension(x, extension)
%
% Default extension: 'wav'
%
% Programmed by Alejandro Osses, ExpORL, 2014
% Last update: 13/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    extension = 'wav';
end

try
    N = length(x);
    M = length(extension);
    
    if strcmp(x(N-M+1:N),extension);
        y = x(1:N-M-1); 
    else
        y = x;
    end
catch
    warning('String seems to be very short, probably it is without extension');
    y = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
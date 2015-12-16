function y = meanrms(x,fs);
%   function y = meanrms(x,fs);
% Calculates the mean rms of a signal using a window of 50 msec

% windowlength = 0.05*fs;
% 
% t = zeros(windowlength,floor(length(x)/windowlength));
% t(:) = x(1:windowlength*floor(length(x)/windowlength));
% 
% y = mean(sqrt(mean(t.^2)));

if (nargin>1)
    warning('Meanrms used');
end

y = sqrt(mean(x.^2));
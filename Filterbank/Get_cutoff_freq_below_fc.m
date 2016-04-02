function f_3dB = Get_cutoff_freq_below_fc(f,h,fc)
% function f_3dB = Get_cutoff_freq_below_fc(f,h,fc)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: r20150814_update.m, r20160316_update_Antoine.m
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 14/08/2015
% Last update on: 14/08/2015 
% Last use on   : 14/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    fs = 44100;
    [b, a] = butter(5, 0.6,'high');
    % Determine frequency response
    [h, w] = freqz(b, a, 2048);
    f = w/pi*(fs/2);
end

if nargin == 3
    idx = max( find( f<= fc ) );
    h = h(1:idx);
    f = f(1:idx);
end

% linear approximation of 3 dB cutoff frequency
ind = find(abs(h) > sqrt(1/2)*max(abs(h)), 1, 'first');
%slope = (abs(h(ind)) - abs(h(ind - 1))) / (f(ind) - f(ind - 1));
f_3dB = f(ind); % ( sqrt(1/2) - abs(h(ind - 1)) + slope * f(ind - 1) ) / slope;
% check result

if nargout == 0
    figure; 
    plot(f,abs(h))
    hold on;
    plot(f_3dB, sqrt(1/2), 'rx');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

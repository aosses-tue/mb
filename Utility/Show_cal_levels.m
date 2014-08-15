function y = Show_cal_levels
% function y = Show_cal_levels
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 14/08/2014 % Update this date manually
% Last use on   : 14/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

callevels = {   'AM tb    ',    '0 dBFS',   '100 dB SPL'; ...
                'fastl2007',    '0 dBFS',   ' 90 dB SPL'};

for i=1:size(callevels,1)
    fprintf('source: %s, calibration level: %s = %s \n',callevels{i,1},callevels{i,2},callevels{i,3})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [h stats] = Boxplot(x)
% function [h stats] = Boxplot(x)
%
% 1. Description:
%       Same as boxplot, but returning all plotted parametes in the struct
%       'stats'
% 
%       h - boxplot handles
%       stats - struct containing: Group, Median, Percentile75th, Percentile25th
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 19/06/2014
% Last update: 19/06/2014 % Update this date manually
% Last used: 23/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = boxplot(x);
grid on;

% get(H(:,3),'tag')

try
    htmp = findobj(h,'tag','Outliers');
    xdata = get(htmp,'XData');
    ydata = get(htmp,'YData');
end

try
    htmp = findobj(h,'tag','Median');
    % stats.Group     = reverse_matrix( mean( cell2mat( get(htmp,'XData') )' ) );
    % stats.Median    = reverse_matrix( mean( cell2mat( get(htmp,'YData') )' ) ); % same value 2 times
    stats.Group     = mean( cell2mat( get(htmp,'XData') )' );
    stats.Median    = mean( cell2mat( get(htmp,'YData') )' ); % same value 2 times
end

try
    htmp = findobj(h,'tag','Upper Whisker');
    stats.Percentile75th = mean( cell2mat( get(htmp,'YData') )' ); % same value 2 times
end

try
    htmp = findobj(h,'tag','Lower Whisker');
    stats.Percentile25th = mean( cell2mat( get(htmp,'YData') )' ); % same value 2 times
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
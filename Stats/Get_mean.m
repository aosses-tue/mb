function [m,s,ci] = Get_mean(data)
% function [m,s,ci] = Get_mean(x)
%
% Get means 'm', standard deviation 's' and confidence interval 'ci' 
% ignoring NaN data.
%
% Programmed by Matthias Milczynski, ExpORL, KU Leuven, Belgium 2008-2011
% Edited by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on : 18/08/2014
% Last use on    : 18/08/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [m,s,ci] = calcMeanStd(data)

m = nan(1, size(data, 2)); % zeros
s = nan(1, size(data, 2)); % zeros
ci = nan(1, size(data, 2)); %zeros
for j=1:size(data, 2)
    d = data(:, j);
    idx = find(isnan(d)); %-1;
    if length(idx > 0)
        warning(sprintf('Excluding from mean and std calculation %.0f NaN data point(s)',length(idx)));
    end
    d(idx) = [];
    if ~isempty(d)
        %if obj.interval
        m(1, j) = mean(d);
        s(1, j) = std(d);
        ci(1, j) = confidenceInterval(d);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

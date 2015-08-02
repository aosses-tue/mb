function [m,s,ci,iq] = Get_mean(data)
% function [m,s,ci,iq] = Get_mean(data)
%
% 1. Description:
%       Gets the typical statistics ignoring NaN data:
%           'm' = mean
%           's' = standard deviation 
%           'ci' = confidence interval
%           'iq' = Interquartile range
%
% Programmed by Matthias Milczynski, ExpORL, KU Leuven, Belgium 2008-2011
% Edited by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on : 02/08/2015
% Last use on    : 02/08/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [m,s,ci] = calcMeanStd(data)

m = nan(1, size(data, 2)); % zeros
s = nan(1, size(data, 2)); % zeros
ci = nan(1, size(data, 2)); %zeros
iq = nan(3, size(data, 2)); %zeros

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
        percL = 25;
        percU = 75;
        [Median,errorL,errorU] = Prepare_errorbar_perc(d,percL,percU);
        iq(1, j) = Median;
        iq(2, j) = errorL;
        iq(3, j) = errorU;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [m,s,ci] = Get_mean(data)
% function [m,s,ci] = Get_mean(x)
%
% Get means 'm', standard deviation 's' and confidence interval 'ci' 
% ignoring NaN data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [m,s,ci] = calcMeanStd(data)

m = nan(1, size(data, 2)); % zeros
s = nan(1, size(data, 2)); % zeros
ci = nan(1, size(data, 2)); %zeros
for j=1:size(data, 2)
    d = data(:, j);
    idx = find(isnan(d)); %-1;
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
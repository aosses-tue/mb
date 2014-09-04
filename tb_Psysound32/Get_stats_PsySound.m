function [stats Nparams bTimeSeries] = Get_stats_PsySound(tmpObj)
% function [stats Nparams bTimeSeries] = Get_stats_PsySound(tmpObj)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 04/09/2014
% Last update on: 04/09/2014 % Update this date manually
% Last use on   : 04/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = [];
Nparams = length(tmpObj);
stats.min   = nan(1,Nparams);
stats.max   = nan(1,Nparams);
stats.mean  = nan(1,Nparams);
stats.percentiles = nan(99,Nparams);

for i = 1:Nparams
    tmpstats = get(tmpObj{1,i},'Stats');
    
    if length(tmpstats) ~= 0
        bTimeSeries = 1;
        stats.min(i) = get(tmpstats,'min');
        stats.max(i) = get(tmpstats,'max');
        stats.mean(i) = get(tmpstats,'mean');
        stats.percentiles(:,i) = get(tmpstats,'percentiles');
    else
        bTimeSeries = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

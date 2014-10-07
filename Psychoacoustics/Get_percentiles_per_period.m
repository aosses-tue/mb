function pXX = Get_percentiles_per_period(x,samples_per_period)
% function y = Get_percentiles_per_period(x,samples_per_period)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 24/09/2014
% Last update on: 24/09/2014 % Update this date manually
% Last use on   : 24/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = buffer(x,samples_per_period,0);
y = y(:,1:end-1); % Delete last period (probably it is truncated)

for i = 5:5:95
    str_exp = sprintf('pXX.p%.0f = percentile(y,%.0f);',i,i);
    eval(str_exp);
end

pXX.timeseries = y;

if nargout == 0
    fprintf('Computed percentiles (mean value), considering %.0f periods (p5, 50 and 95):\n',length(pXX.p5));
    disp(Num2str(mean(pXX.p5)))
    disp(Num2str(mean(pXX.p50)))
    disp(Num2str(mean(pXX.p95)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

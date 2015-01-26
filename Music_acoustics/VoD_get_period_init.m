function params = VoD_get_period_init(bOption)
% function params = VoD_get_period_init(bOption)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%   Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 28/07/2014
% Last update on: 28/07/2014 % Update this date manually
% Last used on  : 28/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bOption = 0; % current values
end

switch bOption
    case 0
        
        params.ti_measured = [0.2616 0.0927 0.0494 0.0961]; 
        params.ti_model    = [0.6023 0.3934 0.2961 0.2688]; % period starts when model starts, but one period later (to avoid problem with offset due to HPF)
        
    case 1 % up to 25/07/2014
        
        % Columns represent different takes, rows different modes. At the specified
        % times 'ti' in seconds the first 'zero' values were found. These values were
        % obtained using VoD_get_initial_time_measured followed by a manual time 
        % correction (in the respective figure, 'Data cursor' was used at the 
        % specified 'stats.tmin' times):
        params.ti = [ NaN     0.1122   NaN     NaN; ...   % take 2: nicer
                      NaN     NaN      0.0273  NaN; ...   % take 3: nicer
                      NaN     NaN      NaN     0.0131; ...% take 4: nicer
                      NaN     NaN      0.1727  NaN];      % take 3: nicer
        for i = 1:4
            params.ti_measured(i) = params.ti(i,params.last_take(i));
        end
        
        params.ti_model      = [0.4796 0.1178 0.2545 0.0732]; 
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function r20160404_update_proc_ICRA_APEX
% function r20160404_update_proc_ICRA_APEX
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 01/04/2016
% Last update on: 01/04/2016 
% Last use on   : 01/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_where = ['D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-5-Multiprocedure-auto\S00-Initial-test-AO' delim];

f = {'piano_multi_result-C2-test-1.apr', 'piano_multi_result-A4-test-1.apr', 'piano_multi_Csh5-test-1.apr'};

h = [];
for i = 1:length(f)
    
    quick_staircases([dir_where f{i}]);
    h(end+1) = Figure2paperfigureT(gcf,2,2);
    
    Saveas(h(end),['staircase-' num2str(i)]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

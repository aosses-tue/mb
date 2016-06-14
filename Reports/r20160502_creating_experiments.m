function y = r20160502_creating_experiments
% function y = r20160502_creating_experiments
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 02/05/2016
% Last update on: 02/05/2016 
% Last use on   : 02/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
bExperiments_op_0502 = 1;

if bExperiments_op_0502;
    opts.dirstimuli = 'Stimuli-20150502';
    Experiments = [17.1 17.3 17.5 17.6];
    do_skip = 0;
    bAdaptive = 0;
    Create_Piano_ICRA_multi_20160418(Experiments,do_skip,bAdaptive,opts);
    % 'PILOT-CONST-AB'

    Experiments = [17   17.2 17.4 17.7];
    do_skip = 0;
    bAdaptive = 0;
    Create_Piano_ICRA_multi_20160418(Experiments,do_skip,bAdaptive,opts);
    % 'PILOT-CONST-BA'

    Experiments = [     17.3 17.5 17.6];
    do_skip = 1;
    bAdaptive = 1;
    Create_Piano_ICRA_multi_20160418(Experiments,do_skip,bAdaptive,opts);
    % 'PILOT-ADAPT-AB'

    Experiments = [     17.2 17.4 17.7];
    do_skip = 1;
    bAdaptive = 1;
    Create_Piano_ICRA_multi_20160418(Experiments,do_skip,bAdaptive,opts);
    % 'PILOT-ADAPT-BA'
end


opts.dirstimuli = 'Stimuli-20150503';
Experiments = [17.1 17.3 17.5 17.6];
do_skip = 0;
bAdaptive = 0;
Create_Piano_ICRA_multi_20160418(Experiments,do_skip,bAdaptive,opts);
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

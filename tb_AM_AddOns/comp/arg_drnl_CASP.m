function definput=arg_drnl_CASP(definput)
% function definput=arg_drnl_CASP(definput)
%
% Created/edited by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Original filename: arg_drnl.m
% Created on    : 15/08/2015
% Last update on: 15/08/2015 
% Last use on   : 15/08/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % 0. Parameters as in the Gammatone filterbank:
    definput.keyvals.flow     = 80;
    definput.keyvals.fhigh    = 8000;
    definput.keyvals.basef    = [];
    definput.keyvals.bwmul    = 1;

    % parameters according to Lopez-Poveda and Meddis 2001 but updated to Jepsen 2008
    %% Linear part:
    definput.keyvals.lin_ngt  = 2; % Number of Gammatone filters
    definput.keyvals.lin_nlp  = 4; % Number of LP filters
    definput.keyvals.lin_fc   = [-0.06762 1.01679];
    definput.keyvals.lin_bw   = [  .03728  .75   ]; % updated
    definput.keyvals.lin_gain = [ 4.20405 -.47909];
    definput.keyvals.lin_lp_cutoff = [-0.06762 1.01 ]; % updated
  
  %% Non-linear part:
    definput.keyvals.nlin_ngt_before  = 2; % first-order. In paper these were 3-cascade filters
    definput.keyvals.nlin_ngt_after   = 2; % first-order. In paper these were 3-cascade filters
    definput.keyvals.nlin_nlp         = 1; % second-order. In paper these were 3-cascade filters
    definput.keyvals.nlin_fc_before = [-0.05252 1.01650];
    definput.keyvals.nlin_fc_after  = [-0.05252 1.01650];
    definput.keyvals.nlin_bw_before = [-0.03193  .77   ]; % updated
    definput.keyvals.nlin_bw_after  = [-0.03193  .77   ]; % updated
    definput.keyvals.nlin_lp_cutoff = [-0.05252 1.01650];
 
    % broken stick non-linearity frequencies below 1000 Hz:
    definput.keyvals.nlin_a = [1.40298 .81916 ];
    definput.keyvals.nlin_b = [1.61912 -.81867 ];
    definput.keyvals.nlin_c = [log10(.25) 0];
    
    % definput.keyvals.compresslimit = [1500];
   
    %% Other parameters of the CASP model
    % The first in the cell is the default setting:
    definput.flags.outerear = {'outerear'};
    definput.flags.middleear={'jepsenmiddleear','middleear','nomiddleear'};
    definput.flags.path = {'bothparts','linonly','nlinonly'};
    definput.flags.ihctype = {'ihc_jepsen'};
    definput.flags.modfiltertype = {'modfilterbank','lowpass'};
    definput.flags.resample_intrep = {'resample_intrep','noresample_intrep'};

    % This parameter set is not supported anymore, as there is no evince as
    % to whether or not this is actually the dataset that was used in the paper.  
    definput.groups.jepsen2008={...
                                'lin_bw',          [ .03728   .78563 ],... % 'lin_lp_cutoff',   [-0.06762 1.01 ],...   % 'nlin_bw_before' , [-0.03193  .77 ],... %                                 'compresslimit', 1500, ... %                                 'jepsenmiddleear', ...
                                'nlin_ngt_before', 2 };

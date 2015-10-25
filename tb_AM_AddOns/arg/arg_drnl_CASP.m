function definput = arg_drnl_CASP(definput)
% function definput = arg_drnl_CASP(definput)
%
% Created/edited by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Original filename: arg_drnl.m
% Created on    : 15/08/2015
% Last update on: 27/08/2015 
% Last use on   : 27/08/2015
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

% How obtained:
%       'p0' = -3*MM
%       'm' = MM = ( From_dB(2.5)-From_dB(1) )/(log10(4000)-log10(1000))
% definput.keyvals.gain_after_compression = [-1.0539 0.3513]; % this is the closest, but 1k seems still biased
definput.keyvals.gain_after_compression = [0 0]; % this is the closest, but 1k seems still biased

%% Non-linear part:
definput.keyvals.nlin_ngt_before  = 2; % first-order. In paper these were 3-cascade filters
definput.keyvals.nlin_ngt_after   = 2; % first-order. In paper these were 3-cascade filters
definput.keyvals.nlin_nlp         = 1; % second-order. In paper these were 3-cascade filters
definput.keyvals.nlin_fc_before = [-0.05252 1.01650];
definput.keyvals.nlin_fc_after  = [-0.05252 1.01650];
definput.keyvals.nlin_bw_before = [-0.03193  .77   ]; % updated
definput.keyvals.nlin_bw_after  = [-0.03193  .77   ]; % updated
definput.keyvals.nlin_lp_cutoff = [-0.05252 1.01650];

% definput.keyvals.compresslimit = [1500];
% broken stick non-linearity frequencies below compresslimit of 1000 [Hz]:
definput.keyvals.nlin_a = [1.40298  .81916];
definput.keyvals.nlin_b = [1.61912 -.81867];
definput.keyvals.nlin_c = [log10(.25) 0];

% broken stick non-linearity frequencies above compresslimit of 1000 [Hz]:
definput.keyvals.nlin_a_above = [ 4.00471 0];
definput.keyvals.nlin_b_above = [-0.98015 0];
% definput.keyvals.nlin_a_above = [1.40298  .81916];
% definput.keyvals.nlin_b_above = [1.61912 -.81867];

definput.keyvals.maxfreqs = [86.9    257.4   520    1055.8       3982.6  7819.2];
definput.keyvals.maxouts = [-21.5255 -21.8778 -19.2899 -22.0003 -35.3780 -45.1737]; %[-26.0464 -20.8640 -17.9071 -19.3908 -35.5670 -45.7842];
                           
% figure; plot(maxfreqs,maxouts);
    
% Gain after DRNL:
definput.keyvals.gain_after_drnl = 50;

%% Other parameters of the CASP model
% The first in the cell is the default setting:
definput.flags.outerear         = {'outerear','noouterear'};
definput.flags.middleear        = {'middleear','jepsenmiddleear','nomiddleear'};
definput.flags.path             = {'bothparts','linonly','nlinonly'};
definput.flags.absolutethreshold = {'noabsolutethreshold','absolutethreshold'};
definput.flags.ihctype          = {'ihc_jepsen'};
definput.flags.modfiltertype    = {'modfilterbank','lowpass'};
definput.flags.resample_intrep  = {'noresample_intrep','resample_intrep'};

definput.groups.lopezpoveda2001 = {... % taken from Lopez-Poveda, Table III
                            'lin_lp_cutoff'  , [-0.06762 1.01679], ...
                            'lin_bw'         , [ 0.03728  .78563], ...
                            'nlin_bw_before' , [-0.03193  .77426], ...
                            'nlin_bw_after'  , [-0.03193  .77426], ...
                            'nlin_a_above'   , [1.40298  .81916], ... % same as 'nlin_a'
                            'nlin_b_above'   , [1.61912 -.81867], ... % same as 'nlin_b'
                            'nlin_ngt_before', 3, ...   % Fig. 3.a
                            'nlin_ngt_after' , 3, ...   % Fig. 3.a
                            'nlin_nlp'       , 3};      % Fig. 3.a
                           
definput.groups.jepsen2008={...
                            'lin_bw',          [ .03728   .78563 ],... % 'lin_lp_cutoff',   [-0.06762 1.01 ],...   % 'nlin_bw_before' , [-0.03193  .77 ],... %                                 'compresslimit', 1500, ... %                                 'jepsenmiddleear', ...
                            'nlin_ngt_before', 2, ...
<<<<<<< HEAD
                            'gain_after_drnl', 23.86}; % Added by AO
=======
                            'gain_after_drnl', 50};
>>>>>>> 374d1a70ef4d0fc079fce42790cef14d3575841d
    

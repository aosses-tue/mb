function definput=arg_ihcenvelope(definput)
% function definput=arg_ihcenvelope(definput)
%
%  1. Description:
%   'ihc_jepsen added by AO'
%
% Last update on: 28/06/2015 
% Last use on   : 11/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  definput.flags.ihctype={'nodefault','ihc_bernstein','ihc_breebaart','ihc_dau', ...
                          'hilbert','ihc_lindemann','ihc_meddis','ihc_jepsen'};

  definput.keyvals.minlvl=[];
  definput.keyvals.minlvl=[];
  definput.keyvals.ihc_filter_order=5; % added in v. 0.9.7 and only used as parameter
                                       % in flags.do_ihc_breebaart

%   Url: http://amtoolbox.sourceforge.net/doc/comp/arg_ihcenvelope.php

function definput=arg_auditoryfilterbank(definput)
% function definput=arg_auditoryfilterbank(definput)
%
% Last used on: 27/06/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.flow=80;
definput.keyvals.fhigh=8000;
definput.keyvals.basef=[];
definput.keyvals.bwmul=1;

definput.groups.gtf_dau = {'basef',1000};

%   Url: http://amtoolbox.sourceforge.net/doc/comp/arg_auditoryfilterbank.php

% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.7


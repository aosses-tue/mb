function definput=arg_auditoryfilterbank(definput)
% function definput=arg_auditoryfilterbank(definput)
%
% 1. Description:
%       Default parameters for Gammatone filter bank.
% 
%       gtf_dorp added by AO
% 
% Last used on: 27/06/2015
% Last edited on: 01-07-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.flow=80;
definput.keyvals.fhigh=8000;
definput.keyvals.basef=[];
definput.keyvals.bwmul=1;

definput.groups.gtf_dau = {'basef',1000};
definput.groups.gtf_dorp = {'flow',168,'fhigh',1840};  

%   Url: http://amtoolbox.sourceforge.net/doc/comp/arg_auditoryfilterbank.php

% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.7


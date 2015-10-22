function definput=arg_eicell(definput)

  definput.keyvals.tc    = 30e-3;   % Temporal smoothing constant
  definput.keyvals.rc_a  = 0.1;     % Range compression parameter 'a' 
  definput.keyvals.rc_b  = .00002;  % Range compression parameter 'b'
  definput.keyvals.ptau  = 2.2e-3;  % time constant for p(tau) function

%
%   Url: http://amtoolbox.sourceforge.net/doc/comp/arg_eicell.php

% Copyright (C) 2009-2015 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.7
%

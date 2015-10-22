function definput=arg_lindemann1986bincorr(definput)

definput.keyvals.c_s   = 0.3;
definput.keyvals.w_f   = 0.035;
definput.keyvals.M_f   = 6;
definput.keyvals.T_int = 5;
definput.keyvals.N_1   = 1;

definput.groups.stationary={'T_int',Inf,'N_1',17640'};

%
%   Url: http://amtoolbox.sourceforge.net/doc/arg/arg_lindemann1986bincorr.php

% Copyright (C) 2009-2015 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.7
%


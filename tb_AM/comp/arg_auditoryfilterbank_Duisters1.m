function definput=arg_auditoryfilterbank_Duisters1(definput)
% function definput=arg_auditoryfilterbank_Duisters1(definput)
%
% Last used on: 27/06/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flow    = audtofreq( 1,'erb');
fhigh   = audtofreq(40,'erb');
definput.keyvals.flow   = flow;
definput.keyvals.fhigh  = fhigh;
definput.keyvals.basef  = [];
definput.keyvals.bwmul  = 1;

%   Url: http://amtoolbox.sourceforge.net/doc/comp/arg_auditoryfilterbank.php

% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

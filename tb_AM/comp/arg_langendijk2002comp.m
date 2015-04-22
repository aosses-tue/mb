function definput=arg_langendijk2002comp(definput)
% function definput=arg_langendijk2002comp(definput)
% 
%  1. Description:
%       It defines the default input arguments for plotting the experimental
%       data of Langendijk.
% 
%   Url: http://amtoolbox.sourceforge.net/doc/comp/arg_langendijk2002comp.php
% 
% Copyright (C) 2009-2014 Peter L. Sondergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option) 
% any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
% for more details.
%
% You should have received a copy of the GNU General Public License along 
% with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Edited by: Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last edited on: 22/04/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.cp={'std','xcorr'};  

definput.keyvals.s=2;
definput.keyvals.do=0;
  
end

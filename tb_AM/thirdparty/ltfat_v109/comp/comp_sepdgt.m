function [coef]=comp_sepdgt(f,g,a,M)
%COMP_SEPDGT  Separable DGT
%   Usage:  c=comp_sepdgt(f,g,a,M);
%  
%   This is a computational routine. Do not call it directly.
%
%   See help on DGT.
%
%   Url: http://ltfat.sourceforge.net/doc/comp/comp_sepdgt.php

% Copyright (C) 2005-2014 Peter L. Søndergaard <soender@users.sourceforge.net>.
% This file is part of LTFAT version 1.4.4
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

%   AUTHOR : Peter L. Søndergaard.

L=size(f,1);
Lwindow=size(g,1);

if Lwindow<L
   % Do the filter bank algorithm
   coef=comp_dgt_fb(f,g,a,M);
else
   % Do the factorization algorithm
   coef=comp_dgt_long(f,g,a,M);
end;

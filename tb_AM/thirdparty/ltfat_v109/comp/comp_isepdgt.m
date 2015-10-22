function [f]=comp_isepdgt(coef,g,L,a,M)
%COMP_ISEPDGT  Separable IDGT.
%   Usage:  f=comp_isepdgt(c,g,L,a,M);
%       
%   This is a computational routine. Do not call it directly.
%
%   Input must be in the M x N x W format.
%
%   See also: idgt
%
%   Url: http://ltfat.sourceforge.net/doc/comp/comp_isepdgt.php

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
%   TESTING: OK
%   REFERENCE: OK

Lwindow=size(g,1);

    if L==Lwindow
        % Do full-window algorithm.
        % coef=reshape(coef,M,prod(size(coef))/M);
        
        % Get the factorization of the window.
        %gf = comp_wfac(g,a,M);      

        % Call the computational subroutine.
        f  = comp_idgt_long(coef,g,L,a,M);
    else
        %coef=reshape(coef,M,prod(size(coef))/M);
        % Do filter bank algorithm.
        % Call the computational subroutine.

        f=comp_idgt_fb(coef,g,L,a,M);
    end;

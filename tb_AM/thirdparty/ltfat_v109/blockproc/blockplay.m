function blockplay(f)
%BLOCKPLAY Schedules block to be played
%   Usage: blockplay(L)
%       
%   Input parameters:
%      f    : Samples.
%
%   Function schedules samples in f to be played. Since playrec handles
%   playing and recording in a single command, the actual relay of samples
%   to playrec is done in the next call of BLOCKREAD.
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/blockplay.php

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


source = block_interface('getSource');

if strcmp(source,'rec')
   error('%s: Blocks cannot be played in the rec only mode.',upper(mfilename));
end

block_interface('setToPlay',f);


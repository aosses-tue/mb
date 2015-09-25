function [mue corrmue] = optimaldetector(ir_stim,template,fs)
% function [mue corrmue] = optimaldetector(ir_stim,template,fs)
% 
%   1. Description:
%       OPTIMALDETECTOR  Generic optimal detector for the CASP and Breebaart models
%
%       This is a correlation-based optimal detector for a signal known 
%       exactly see Green & Swets (1966) for more details.
%
%       Url: http://amtoolbox.sourceforge.net/doc/modelstages/optimaldetector.php
%
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    normmethod = 2; % default in AMT
else
    normmethod = 1;
end

corrmue = ir_stim.*template;

switch normmethod
    case 1
        optfactor = 1/fs;
        mue = sum(corrmue(:))*optfactor;
    case 2
        optfactor = sqrt(numel(corrmue)); % default in AMT
        % Take mean over all dimensions of internal representation and correct for
        % optimalityfactor.
        mue = mean(corrmue(:))*optfactor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
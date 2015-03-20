function y = rms(insig,options)
%RMS returns the root mean square of the signal
%
%   Usage: y = rms(insig,[options]);
%
%   Input parameters:
%       insig       - signal for which to calculate the RMS value
%       options     - if 'ac' is given, only the AC component of the signal
%                     is used for RMS calculation
%
%   Output parameters:
%       y           - RMS value of insig
%
%   RMS(insig) computes the RMS (Root Mean Square) value of a finite
%   sampled signal sampled at a uniform sampling rate.
%
%   RMS(x,'ac') does the same, but considers only the AC component of the
%   signal (i.e. the mean is removed).
%
%   If the input is a matrix or ND-array, the RMS is computed along the first
%   dimension, and a vector of values is returned.
%
%   The RMS value of a signal x of length N is computed by
%
%                        1  N
%     rms(insig) = sqrt( - sum insig(n)^2 )
%                        N n=1
%
%   see also: db

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************
% Edited by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Edited on     : 25/05/2014
% Last used on  : 01/06/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
try
    narginchk(nargmin,nargmax);
catch
    warning('narginchk not done')
end
isargnumeric(insig);
if (nargin==1) || (~ischar(options))
  options='';
end


%% ===== Computation =====================================================
% It is better to use 'norm' instead of explicitly summing the squares, as
% norm (hopefully) attempts to avoid numerical overflow.
y = zeros(1,size(insig,2));
switch(lower(options))
    case 'ac'
        for ii=1:size(insig,2)
            y(1,ii) = norm(insig(:,ii)-mean(insig(:,ii)))/sqrt(length(insig(:,ii)));
        end
    otherwise
        for ii=1:size(insig,2)
            y(1,ii) = norm(insig(:,ii))/sqrt(length(insig(:,ii)));
        end
end

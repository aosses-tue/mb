% DPRIME  - Signal-detection theory sensitivity measure.
%
%  d = dprime(pHit,pFA)
%  [d,beta] = dprime(pHit,pFA)
%
%  PHIT and PFA are numerical arrays of the same shape.
%  PHIT is the proportion of "Hits":        P(Yes|Signal)
%  PFA is the proportion of "False Alarms": P(Yes|Noise)
%  All numbers involved must be between 0 and 1.
%  The function calculates the d-prime measure for each <H,FA> pair.
%  The criterion value BETA can also be requested.
%  Requires MATLAB's Statistical Toolbox.
%
%  References:
%  * Green, D. M. & Swets, J. A. (1974). Signal Detection Theory and
%    Psychophysics (2nd Ed.). Huntington, NY: Robert Krieger Publ.Co.
%  * Macmillan, Neil A. & Creelman, C. Douglas (2005). Detection Theory:
%    A User's Guide (2nd Ed.). Lawrence Erlbaum Associates.
%  
%  See also NORMINV, NORMPDF.

% Original coding by Alexander Petrov, Ohio State University.
% $Revision: 1.2 $  $Date: 2009-02-09 10:49:29 $
%
% Part of the utils toolbox version 1.1 for MATLAB version 5 and up.
% http://alexpetrov.com/softw/utils/
% Copyright (c) Alexander Petrov 1999-2006, http://alexpetrov.com
% Please read the LICENSE and NO WARRANTY statement:

% GNU Public License for the UTILS Toolbox
%  
% Alex Petrov, December 2006
% 
% Downloaded on: 08/10/2014 
% Downloaded from: https://sccn.ucsd.edu/svn/software/eeglab/functions/miscfunc/dprime.m
% ==============================================================================

function [d,beta] = dprime(pHit,pFA)

%-- Convert to Z scores, no error checking
zHit = norminv(pHit) ;
zFA  = norminv(pFA) ;

%-- Calculate d-prime
d = zHit - zFA ;

%-- If requested, calculate BETA
if (nargout > 1)
  yHit = normpdf(zHit) ;
  yFA  = normpdf(zFA) ;
  beta = yHit ./ yFA ;
end

%%  Return DPRIME and possibly BETA
%%%%%% End of file DPRIME.M

function theRMS = calcRMS(input,dim)
%calcRMS   Calculate the root mean squared value.
% 
%USAGE
%   theRMS = calcRMS(input)
%   theRMS = calcRMS(input,dim)
% 
%INPUT ARGUMENTS
%    input : N-dimensional input signal 
%      dim : dimension along which the RMS should be calculated
% 
%OUTPUT ARGUMENTS
%   theRMS : (N-1)-dimensional RMS
% 
%   See also calcRSS and setLevelRMS.

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to
%   
%   Author  :  Tobias May, © 2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/10/15
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% CALCULATE THE ROOT MEAN SQAURED VALUE
% 
% 
% Calculate the RMS
if nargin < 2 || isempty(dim); 
    theRMS = sqrt(mean(input .* conj(input)));
else
    theRMS = sqrt(mean(input .* conj(input), dim));
end

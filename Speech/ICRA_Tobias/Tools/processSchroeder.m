function out = processSchroeder(input)
%processSchroeder   Flip the sign of 50 per cent of the input samples.
%   This process will maintain the modulation properties of the input, but
%   renders the signal (e.g. speech) unintelligible. The output signal has
%   a flat (white) spectrum [1,2].   
% 
%USAGE
%      out = processSchroeder(input)
%
%INPUT ARGUMENTS
%    input : input signal [nSamples x nChannels x ... ]
% 
%OUTPUT ARGUMENTS
%   output : output signal [nSamples x nChannels x ... ]
% 
%   See also randPerturbation and randomizePhase.
% 
%REFERENCES
%   [1] M. R. Schroeder, "Reference Signal for Signal Quality Studies", The
%       Journal of the Acoustical Society of America, 44(6), pp. 1735-1736,
%       1968.    
%  
%   [2] W. A. Dreschler, H. Verschuure, C. Ludvigsen and S. Westermann,
%       "ICRA Noises: Artificial Noise Signals with Speech-like Spectral
%       and Temporal Properties for Hearing Instrument Assessment",
%       International Journal of Audiology, 40(3), pp.148-157, 2001. 


%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, (c) 2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/02/09
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 1
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% PERFORM SCHROEDER PROCESSING
% 
% 
% Determine size of input signal
dim = size(input);

% Randomly change the sign of 50 % of the samples
mod = 2 * (rand(dim) > 0.5) - 1;

% Flip sign
out = input .* mod;

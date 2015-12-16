function out = frameData(input,blockSize,stepSize,win,bZeroPad)
%frameData   Segment input signal into overlapping frames.
% 
%USAGE
%      frames = frameData(input,blockSize,stepSize,win,bZeroPad);
% 
%INPUT ARGUMENTS
%       input : input signal [nSamples x 1]
%   blockSize : frame size in samples
%    stepSize : step size between frames in samples 
%               (default, stepSize = round(0.5 * blockSize))
%         win : string or vector defining the analysis window 
%               (default, win = 'rectwin')
%    bZeroPad : zero-pad input signal to fill last frame 
%               (default, bZeroPad = true)
% 
%OUTPUT ARGUMENTS
%      frames : frame-based signal [blockSize x nFrames]
%
%   See also frames2audio.
% 
%EXAMPLE
%   % Use frames2audio to reconstruct a time domain signal from a matrix of
%   % overlapping frames:
% 
%   % Input signal 
%   input = randn(1E3,1); 
% 
%   % Select a square-root and periodic analysis window
%   winAnalysis = sqrt(hamming(256,'periodic'));
% 
%   % Segment input into frames of 256 samples with 128 samples overlap
%   frames = frameData(input,256,128,winAnalysis,true);
% 
%   % Reconstruct signal and visualize synthesis
%   out = frames2audio(frames,128,winAnalysis,'synthesis',true);
%
%   % Trim output signal
%   out = out(1:numel(input));
% 
%   % Calculate MSE error in dB
%   calcMSE(input,out)


%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :
%   v.0.1   2009/11/23
%   v.0.2   2010/05/18 automatically determine the number of frames
%   v.0.3   2013/08/19 added zero-padding
%   v.0.4   2015/02/20 added documentation
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(stepSize); stepSize = round(0.5 * blockSize); end
if nargin < 4 || isempty(win);      win      = 'hamming';              end
if nargin < 5 || isempty(bZeroPad); bZeroPad = true;                   end
    
% Determine size of input signal
[nSamples,nChannels] = size(input);

% Check for proper size
if nChannels > 1     
    error('Monaural input is required!')
end


%% ****************************  ZERO-PADDING  ****************************
% 
% 
% Compute number of frames
nFrames = (nSamples-(blockSize-stepSize))/(stepSize);

% Append zeros
if bZeroPad
    % Number of frames (at least one)
    nFrames = ceil(max(1,nFrames));
    
    % Compute number of required zeros
    nZeros = (stepSize * nFrames + (blockSize-stepSize)) - nSamples;
    
    % Pad zeros
    input = [input; zeros(nZeros,1)];
else
    % No zero padding, exclude elementes that do not fit into last frame
    nFrames = max(0,floor(nFrames));
end
 
% Check if window is a window in samples or a string | function handle
if ischar(win) || isa(win,'function_handle')
    win = window(win,blockSize);
else
    if blockSize ~= length(win)    
       error('Mismatch between blockSize and window size.') 
    end
end
        
    
%% ***********************  FRAME-BASED PROCESSING  ***********************
% 
% 
% Framing (MEX processing)
out = frameDataMEX(input,blockSize,stepSize,win,nFrames);


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************
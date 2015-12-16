function out = frames2audio(frames,hopSize,win,scaling,bPlot)
%frames2audio   Reconstruct a continuous signal from overlapping frames.
% 
%USAGE
%       out = frames2audio(frames,hopSize,win,scaling,bPlot)
% 
%INPUT ARGUMENTS
%    frames : 2D matrix with overlapping frames [blockSize x nFrames]
%   hopSize : step size between frames in samples 
%       win : string or vector defining the analysis window 
%   scaling : string specifying scaling of output signal
%                  'none' - no scaling
% 
%              'analysis' - output is scaled by the sum of the analysis
%                           window, assuming that the overlapping frames
%                           have been windowed during the analysis stage.
%                           No additional windowing is performed during the
%                           synthesis stage  
% 
%             'synthesis' - ouput is scaled by the squared sum of the
%                           analysis window, assuming that the frames have
%                           been windowed during the analysis stage. The
%                           supplied window is applied again during the
%                           synthesis stage. If this option is chosen, it
%                           is common to apply a square-root window during
%                           analysis and synthesis, which will result in a
%                           constant scaling factor. In case the
%                           overlapping frames have been modified prior to
%                           reconstruction, the additional synthesis window
%                           can help to reduce edge discontinuities.  
% 
%     bPlot : visualize the synthesis stage (default, bPlot = false)
% 
%   frames2audio(...)  plots the synthesis stage in a new figure.
% 
%OUTPUT ARGUMENTS
%       out : reconstructed signal [nSamples x 1]
%
%   See also frameData.
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
% 
%ACKNOWLEDGEMENT
%   The different scaling options were adopten from frames2vec.m, a
%   function programmed by Kamil Wojcicki. 


%   Developed with Matlab 8.1.0.604 (R2013a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2007-2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2013/08/19 re-written
%   v.0.2   2015/02/20 added documentation
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 4 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 5 || isempty(bPlot); bPlot = false; end

% Determine size of input signal
[blockSize,nFrames,nChannels] = size(frames);

% Check for proper size
if nChannels > 1
    error('Monaural input is required!')
end


%% APPLY SYNTHESIS WINDOW
% 
% 
% Check if win is a window in samples or a string
if ischar(win)
    win = window(win,blockSize);
else
    if blockSize ~= length(win)
        error('Mismatch between blockSize and window size.')
    end
end

% Select scaling method
switch(lower(scaling))
    case {false 'no' 'none'}
        % No scaling
        wRescale  = 0;
        bScale    = false;
        bApplyWin = false;
    case 'analysis'
        % Output is scaled by the sum of the window 
        wRescale  = win;
        bScale    = true;
        bApplyWin = false;
    case 'synthesis'
        % Apply a synthesis window, therefore the output is scaled by the
        % squared sum of the window 
        wRescale  = win.^2;
        bScale    = true;
        bApplyWin = true;
    otherwise
        error(['Reconstruction method ''',lower(scaling),''' is not supported!'])
end

% Apply synthesis window
if bApplyWin
    % Matrix-based windowing
    frames = win(:,ones(1,nFrames)) .* frames;
end


%% PERFORM RECONSTRUCTION
% 
% 
% Determine number of samples
nSamples = blockSize + (nFrames-1) * hopSize;

% Allocate memory 
out = zeros(nSamples,nChannels);

% Block indices
blockIdx = (1:blockSize)';

% Overlap-and-add synthesis
if bScale
    % Allocate memory
   reScale = zeros(nSamples,nChannels);
      
   % Loop over number of frames
   for ii = 1 : nFrames
       % Reconstruct signal by overlap-add-procedure
       out(blockIdx,:,:) = out(blockIdx,:,:) + frames(:,ii,:);
       
       % Accumulate window function
       reScale(blockIdx,:) = reScale(blockIdx,:) + wRescale;
       
       % Set new block limits
       blockIdx = blockIdx + hopSize;
   end
   
   % Prevent division-by-zero
   reScale(reScale==0) = 1;
   
   % Check scaling
   if reScale < 1E-3
       % Consider to incorporate a floor, or select a different window
       % function
       warning('The scaling factor is very small! Consider incorporating a floor, or select a different window function.')
   end
   
   % Scale audio signal according to the accumulated window function
   out = out ./ reScale;
else
    % No scaling
    reScale = 1;
    
    % Loop over number of frames
    for ii = 1 : nFrames
        % Reconstruct signal by overlap-add-procedure
        out(blockIdx,:,:) = out(blockIdx,:,:) + frames(:,ii,:);
        
        % Set new block limits
        blockIdx = blockIdx + hopSize;
    end
end

            
%% PLOT RESULTS
% 
% 
% Open new figure
if bPlot || nargout == 0
    h = figure;
    
    % Block indices
    blockIdx = (1:blockSize)';

    % Loop over number of frames
    for ii = 1 : nFrames
        figure(h);hold on;
        hSig = plot(blockIdx,frames(:,ii,:));
        set(hSig,'color',[0.75 0.75 0.75]);
        hLegend = plot(blockIdx,win,'b','linewidth',1.25);
        
        % Set new block limits
        blockIdx = blockIdx + hopSize;
    end
    figure(h);hold on;grid on;
    if bScale
        hLeg = plot(reScale,'k','LineWidth',2);
        legend([hSig hLegend hLeg],{'signal' 'windows' 'total scaling'})
    else
        legend([hSig hLegend],{'signal' 'windows'})
    end
    xlabel('Number of samples')
    ylabel('Amplitude')
    xlim([1 nSamples])    
end
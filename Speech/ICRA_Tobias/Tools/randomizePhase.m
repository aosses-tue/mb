function out = randomizePhase(input,fsHz,blockSec,stepSec,bShuffle)
%randomizePhase   Randomize phase of input signal in the FFT domain.
%   The processing is performed using an overlap-and-add framework, where
%   the input signal is segmented into overlapping frames. The amount of
%   smoothing can be controlled by adjusting the step size between adjacent
%   frames. An overlap of 87.5 % has been suggested in [1].  
% 
%USAGE
%         out = randomizePhase(input,fsHz)
%         out = randomizePhase(input,fsHz,blockSec,stepSec,bShuffle)
%
%INPUT ARGUMENTS
%       input : input signal [nSamples x 1]
%        fsHz : sampling frequency in Hz
%    blockSec : block size in seconds (default, blockSec = 20E-3)
%     stepSec : step size in seconds  (default, stepSec = blockSec * (1/8))
%    bShuffle : if true, the original phase will be used and shuffled
%               across FFT bins for each frame, rather than replacing the
%               original phase by random values (default, bShuffle = false)
% 
%OUTPUT ARGUMENTS
%      output : output signal [nSamples x 1]
% 
%   See also randPerturbation and processSchroeder.
% 
%REFERENCES
%   [1] W. A. Dreschler, H. Verschuure, C. Ludvigsen and S. Westermann,
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
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(blockSec); blockSec = 20E-3;            end
if nargin < 4 || isempty(stepSec);  stepSec  = blockSec * (1/8); end
if nargin < 5 || isempty(bShuffle); bShuffle = false;            end

% Determine size of input signal
[nSamples,nChannels] = size(input);

% Check for proper size
if nChannels > 1     
    error('Monaural input is required!')
end


%% INITIALIZE PARAMETERS
% 
% 
% Block processing parameters
blockSize = 2 * round(fsHz * blockSec / 2);
stepSize  = round(fsHz * stepSec);
nfft      = pow2(nextpow2(blockSize));
% Ensures perfect reconstruction, window is applied twice (analysis &
% synthesis), therefore a sqrt-window is selected here
win       = sqrt(hanning(blockSize));
bZeroPad  = true;
scaling   = 'synthesis';


%% GO TO FFT DOMAIN
% 
% 
% Frame data via overlap-and-add
frames = frameData(input,blockSize,stepSize,win,bZeroPad);

% FFT domain
spec = time2freq(frames,nfft);

% Determine size of spectrum
[nRealFreqs,nFrames] = size(spec);


%% RANDOMIZE PHASE
% 
% 
if bShuffle
    % Shuffle original phase values
    if rem(nfft,2)
        % Random indices
        [~,rIdx] = sort(rand(nRealFreqs-1,nFrames),1);
        
        % Single-sided power spectrum is odd, only DC is unique
        phase_Mod = angle(spec(rIdx + 1));    
        phase_DC  = angle(spec(1,:));    
        phase_NY  = [];
    else
        % Random indices
        [~,rIdx] = sort(rand(nRealFreqs-2,nFrames),1);
        
        % Single-sided power spectrum is even, DC and Nyquist are unique
        phase_Mod = angle(spec(rIdx + 1));
        phase_DC  = angle(spec(1,:));    
        phase_NY  = angle(spec(end,:));    
    end
else
    % Replace original phase with randomized phase
    if rem(nfft,2)
        % Single-sided power spectrum is odd, only DC is unique
        phase_Mod = 2 * pi * rand(nRealFreqs-1,nFrames) - pi;
        phase_DC  = (rand(1,nFrames) > 0.5) * pi; % Can be either 0 or pi
        phase_NY  = [];
    else
        % Single-sided power spectrum is even, DC and Nyquist are unique
        phase_Mod = 2 * pi * rand(nRealFreqs-2,nFrames) - pi;
        phase_DC  = (rand(1,nFrames) > 0.5) * pi; % Can be either 0 or pi
        phase_NY  = (rand(1,nFrames) > 0.5) * pi; % Can be either 0 or pi
    end
end

% Add spectra, DC and nyquist
phase_Mod = cat(1,phase_DC,phase_Mod,phase_NY);

% Reconstruct spectrum using modified phase
specMod = abs(spec) .* exp(1j * phase_Mod);


%% GO BACK TO TIME DOMAIN
% 
% 
% Go back to time domain (ensure real-valued signal)
time = real(freq2time(specMod,nfft));

% Remove zero-padding introduced by nfft 
time = time(1:blockSize,:);

% Reconstruct time-domain signal via overlap-and-add
out = frames2audio(time,stepSize,win,scaling);

% Trim zero-padding
out = out(1:nSamples);

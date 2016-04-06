function [stft, f, t, info] = stft(x, fs, nfft, wlen, overlap, nwtype,info)
% function: [stft, t, f, info] = stft(x, fs, nfft, wlen, overlap, nwtype,info)
% 
% Short-Time Fourier Transform with MATLAB Implementation.
% Temporal window used: Hamming.
% 
% Input parameters:
%   x    - signal in the time domain, assumes mono signal
%   wlen - length of the hamming window [samples]
%   overlap - percentage of overlapping among segments [%]. Determines hop-size
%   nfft - number of FFT points
%   fs   - sampling frequency, Hz
% 
% Output parameters:
%   f    - frequency vector, Hz
%   t    - time vector, s
%   stft - STFT matrix (time across columns, freq across rows)
% 
% % 3.1 Stand-alone example:
%       stft_example;
% 
% % 3.2 Stand-alone example with default parameters:
%       [x fs] = wavread('stft_audiofile.wav');
%       stft(x,fs);
% 
%       [x fs] = wavread('stft_audiofile.wav');
%       nfft = 8192;
%       wlen = nfft/2;
%       overlap = 87.5; % percent
%       nwtype = 4; % 4 = hamming
%       stft(x,fs, nfft, wlen, overlap, nwtype);
% 
% Author: M.Sc. Eng. Hristo Zhivomirov
% Created on   : 21/12/2013
% Downloaded on: 23/06/2014
% Last used on : 05/04/2016
% Edited by Alejandro Osses, HTI, TU/e, the Netherlands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default parameters:

if nargout == 0
    
    if nargin < 7
        info = [];
    end
    info = Ensure_field(info,'normalise_time_factor',1);
    
end

info = ef(info,'bColourBar',1);

if nargin < 3
    nfft = 4096*2; % number of fft points (recomended to be power of 2)
end

if nargin < 6
    nwtype = 4; % Hamming window
end

if nargin < 4
    wlen = nfft/2; % window length (recomended to be power of 2)
end

if nargin < 5
    overlap = (1/8)*100; % typical values 50, 25, 12.5, 6.25
end

h = wlen - floor(wlen*overlap/100); % hop size (recomended to be power of 2)

if nargin < 2
    fs  = 44100;
    warning('Assuming default value for fs, this might be incorrect');
    pause(2)
end

framedur    = (wlen/fs)*1000; % in mili-seconds
hopdur      = (   h/fs)*1000; % in mili-seconds
disp([mfilename '.m: '])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% represent x as column-vector if it is not
if size(x,2) > 1 & ( size(x,2) > size(x,1) )
    x = x';
end

% length of the signal
xlen = length(x);

% form a periodic hamming window
[win wtype] = Get_window(nwtype,wlen);

info_str = sprintf('Frame size: %.1f [ms], Hop size: %.1f [ms], \noverlap = %.1f[%%], ',framedur,hopdur,overlap);
info_str = [info_str ' window = ' wtype ' ' sprintf(' (fs = %.0f [kHz])',fs/1000)];

disp(info_str);
fprintf('\nIf you want larger or shorter sections modify FFT-size');

% form the stft matrix
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((xlen-wlen)/h);        % calculate the total number of columns
stft = zeros(rown, coln);           % form the stft matrix

if size(x,2) == 2
    stftR = zeros(rown, coln);
end

% initialize the indexes
indx = 0;
col = 1;

% perform STFT
while indx + wlen <= xlen
    % windowing
    xw = x(indx+1:indx+wlen,:).*repmat(win, 1, size(x,2));
    
    % FFT
    X = fft(xw, nfft);
    
    % update the stft matrix
    stft(:,col) = X(1:(rown),1); % only left side
    
    if size(xw,2)==2
        stftR(:,col) = X(1:(rown),2);
    end
    
    % update the indexes
    indx = indx + h;
    col = col + 1;
end

% calculate the time and frequency vectors
t = (wlen/2:h:xlen-wlen/2-1)/fs;
f = (0:rown-1)*fs/nfft;

if nargout >= 4
    info.nfft = nfft;
    info.wlen = wlen;
    info.wtype = wtype;
    info.overlap = overlap;
    info.hope = h;
    info.framedur = framedur;
    info.hopdur = hopdur;
end

if nargout == 0
    
    if exist('stftR','var')
        disp('Stereo signal processed, but STFT plot with only left channel');
    end
    % See stft_example for explanations about the following lines...
    K       = sum(win)/wlen;
    stft    = abs(stft)/wlen/K;
    
    if rem(nfft, 2) % odd nfft excludes Nyquist point
        stft(2:end, :)      = stft(2:end, :).*2;
    else            % even nfft includes Nyquist point
        stft(2:end-1, :)    = stft(2:end-1, :).*2;
    end

    stft = 20*log10(stft + 1e-6);
    
    % figure;
    imagesc(t/info.normalise_time_factor, f, stft);
    colormap('Gray')
    set(gca,'YDir','normal')
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
    
    
    
    if info.normalise_time_factor == 1
        info = Ensure_field(info,'XLabel',['Time [s], ' info_str]);
    else
        info = Ensure_field(info,'XLabel',[sprintf('Time normalised to %.3f [s], ',info.normalise_time_factor) info_str]);
    end
    xlabel(info.XLabel)
    ylabel(sprintf('Frequency [Hz]'))
    info = Ensure_field(info,'txtTitle', 'Amplitude spectrogram of the signal');
    
    if exist('stftR','var')
        info.txtTitle = [info.txtTitle ' (only L-channel plotted)'];
    end
    
    title(info.txtTitle)
    
    if info.bColourBar
        handl = colorbar;
        set(handl, 'FontName', 'Times New Roman', 'FontSize', 14)
        ylabel(handl, 'Magnitude [dB]')
    end
    
    if isfield(info,'scale_dB')
        YTick = get(handl,'YTick');
        idx = find(YTick+info.scale_dB<0);
        YTick(idx) = [];
        YTickLabel = YTick + info.scale_dB;
        set(handl,'YTick',YTick);
        set(handl,'YTickLabel',YTickLabel);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
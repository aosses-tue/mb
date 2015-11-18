function D = pvoc_stft(x, N, w, H, fs)
% function D = pvoc_stft(x, f, w, h, sr)
%
% 1. Description:
%       Short-time Fourier transform.
% 
%       D = stft2(X, N, W, H, SR)                       
% 
%	Returns some frames of short-term Fourier transform of x.  
%       N - number of points of the FFT (default 256);
%       H - hop size in samples (default W/2);  
%       Data is hann-windowed at W pts (N), or rectangular if W=0, or 
%       with W if it is a vector.
%       Without output arguments, will plot like sgram (fs will get
%       axes right, defaults to 8000).
%	See also 'istft.m'.
% dpwe 1994may05.  Uses built-in 'fft'
% $Header: /home/empire6/dpwe/public_html/resources/matlab/pvoc/RCS/stft.m,v 1.4 2010/08/13 16:03:14 dpwe Exp $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2;  N = 256; end
if nargin < 3;  w = N; end
if nargin < 4;  H = 0; end
if nargin < 5;  fs = 8000; end

% expect x as a row
if size(x,1) > 1
  x = x';
end

s = length(x);

if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,N);
  else
    if rem(w, 2) == 0   % force window to be odd-len
      w = w + 1;
    end
    halflen = (w-1)/2;
    halff   = N/2;   % midpoint of win
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win     = zeros(1, N);
    acthalflen = min(halff, halflen);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
  end
else
    win = w;
end

w = length(win);
% now can set default hop
if H == 0
  H = floor(w/2);
end

c = 1;

% pre-allocate output array
d = zeros((1+N/2),1+fix((s-N)/H));

for b = 0:H:(s-N)
    u = win.*x((b+1):(b+N));
    t = fft(u);
    d(:,c) = t(1:(1+N/2))';
    c = c+1;
end;

% If no output arguments, plot a spectrogram
if nargout == 0
    tt = [0:size(d,2)]*H/fs;
    ff = [0:size(d,1)]*fs/N;
    imagesc(tt,ff,20*log10(abs(d)));
    axis('xy');
    xlabel('time [sec]');
    ylabel('freq [Hz]')
    % leave output variable D undefined
else
    % Otherwise, no plot, but return STFT
    D = d;
end

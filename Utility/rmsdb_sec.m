function [y ty stats] = rmsdb_sec(x,fs,dur_buf_s,bOverlap,tstart)
% function [y ty stats] = rmsdb_sec(x,fs,dur_buf_s,bOverlap,tstart)
%
% 1. Description:
%       x - signal to be buffered in N = 50 sample-section
%       t - x-time stamp in seconds
%       tstart - time is seconds from where the minimum is going to be looked for
%       fs - sampling frequency required to convert tstart from seconds to samples
% 
%       y - row vector containing the rms value (dB) of X at a time resolution
%           determined by ty. For this calculation x is normalised.
%       ty - y-time stamp
%       stats - struct with minimum values
%       stats.tmin - minimum during the first second
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       VoD_one_signal;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/06/2014
% Last update on: 25/06/2014 
% Last use on   : 23/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    tstart     = 0;
    tstart_idx = 1;
else
    tstart_idx = tstart*fs;
end

if nargin < 4
    bOverlap = 1;
end

t    = (1:length(x))/fs;

% x    = x/max(abs(x)); % by-passed on 16/03/2016
N    = round(dur_buf_s * fs);

if bOverlap == 1
    Noverlap = round(N/2);
else
    Noverlap = 0;
end

xsec = buffer(x,N,Noverlap);
x    = xsec;

[r,c]=size(x);
if c == 1 %column vector
    y = 10*log10(x'*x/length(x));
elseif r == 1 % row vector
    y = 10*log10(x*x'/length(x));
else % Generic case:
    y = 10*log10(sum(x.*x)/length(x));
end

if nargout >= 2
    ty = buffer(t,N,0);
    ty = ty(1,:);
end

if nargout >=3
    tstart_idx = round(tstart_idx/N);
    bDuringFirstSecond = 0;
    idx = [];
    dB_step = 5;
    dBi = -70-dB_step;
    dB = dBi;
    n_local_minima = 5;
    while (length(idx) <= n_local_minima | bDuringFirstSecond == 0) % It checks the existence of at least n local minima
        
        if length(idx) > n_local_minima % first condition is fulfilled
            dB_step = 3;
        end
        
        dB = dB + dB_step;
        idx = find( y < dB );
        
        if ty( min(idx) ) < 1
            bDuringFirstSecond = 1;
        end
    end
    stats.idx_ty = idx;
    stats.idx_t  = idx*N;
    stats.dB    = dB;
    stats.min   = ty(idx);
    idx         = idx(idx > tstart_idx); % we exclude all samples below index tstart_idx
    stats.tmin  = min(ty(idx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
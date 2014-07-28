function [y ty stats] = rmsdb_sec(x,t,tstart,fs)
% function [y ty stats] = rmsdb_sec(x,t,tstart,fs)
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
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 24/6/2014
% Last update: 25/6/2014 % Update this date manually
% Last used: 25/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    tstart     = 0;
    tstart_idx = 1;
else
    tstart_idx = tstart*fs;
end

x = x/max(abs(x));
N       = 50;
xsec    = buffer(x,N,0);
x       = xsec;

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
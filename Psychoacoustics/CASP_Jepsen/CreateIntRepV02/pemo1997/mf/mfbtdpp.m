% mfbtdpp.m - post processing for modulation filterbank 'mfbtd.m' (Dau et al. 1997) output.
%				  Gets the real part for centerfrequencies <= 10 Hz and the absolute value
%				  otherwise.
%
% Usage: out = mfbtdpp(in,inf,fs)
%
% in    = input matrix from mfbtd.m.
% inf   = center frequency info vector
% fs	  = sampling rate in Hz
%
% [out1,out2, ...,outn] = output matrix
%
% copyright (c) 1999 Universitaet Oldenburg

function out = mfbtdpp(in,inf,fs)

out=in;

for i=1:length(inf)
    if inf(i) <= 10
        out(:,i) = real(out(:,i));
    else
        out(:,i) = abs(out(:,i));
    end
end

% eof








function y = rmsdb_freq(X,fs,N,ti,tf)
% function y = rmsdb_freq(x,fs,N,ti,tf)
%
% 1. Description:
%       Root-Mean-Square value of x, in dB
%
% 2.1 Example 1:
%   [x, Fs] = wavread('Choice.wav'); 
%   rmsdb(x)
% 
% 2.2 Example 2:
%   y = rmsdb('Choice.wav');
% 
% 2.3 Example 3, rms value between 0.1 and 0.2 seconds:
%   [x, fs] = Wavread('Choice.wav'); 
%   ti = 0.1;
%   tf = 0.2;
%   y = rmsdb(x,fs,ti,tf);
%
% Programmed by ExpORL, KU Leuven, Belgium
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 28/10/2014 % Update this date manually
% Last use on   : 31/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(X)
    try
        X = Wavread(X);
    catch
        error('variable x interpreted as char, but no wav file with such a name was found')
    end
end

if nargin < 5
    Nf = length(X);
else
    Nf = round(tf*fs);
end

if nargin < 4
    Ni = 1;
else
    Ni = ceil(ti*fs + 1e-6); % to make min idx equal to 1
end

if nargin < 3
    N = Nf - Ni;
end

[r,c]=size(X);

if c == 1
    y = 10*log10( X(Ni:Nf)'*X(Ni:Nf)/length(X(Ni:Nf)) );
elseif r == 1
    Nf = c;
    y = 10*log10( X(Ni:Nf)*X(Ni:Nf)'/length(X(Ni:Nf)) );
else % Generic case:
    y = 10*log10( sum(X(Ni:Nf,:).*X(Ni:Nf,:))/length(X(Ni:Nf,:)) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function y = rmsdb(x)
% function y = rmsdb(x)
%
% Root-Mean-Square value of x, in dB
%
% Example 1:
%   [x, Fs] = wavread('Choice.wav'); 
%   rmsdb(x)
% 
% Example 2:
%   y = rmsdb('Choice.wav');
%
% Programmed by ExpORL
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 21/07/2014 % Update this date manually
% Last use on   : 21/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(x)
    try
        x = Wavread(x);
    catch
        error('variable x interpreted as char, but no wav file with such a name was found')
    end
end

[r,c]=size(x);
if c == 1
    y = 10*log10(x'*x/length(x));
elseif r == 1
    y = 10*log10(x*x'/length(x));
else % Generic case:
    y = 10*log10(sum(x.*x)/length(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function y = rmsdbA(x,fs)
% function y = rmsdbA(x,fs)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/08/2014
% Last update on: 21/08/2014 % Update this date manually
% Last use on   : 21/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(x)
    try
        [x fs] = Wavread(x);
    catch
        error('variable x interpreted as char, but no wav file with such a name was found')
    end
end

try
    [b a] = adsgn(fs);
catch
    warning('Assuming fs = 44100 Hz');
    [b a] = adsgn(44100);
end

xA = filter(b,a,x);
[r,c]=size(xA);
if c == 1
    y = 10*log10(xA'*xA/length(xA));
elseif r == 1
    y = 10*log10(xA*xA'/length(xA));
else % Generic case:
    y = 10*log10(sum(xA.*xA)/length(xA));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

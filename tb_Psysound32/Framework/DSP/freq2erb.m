function erb = freq2erb(freq)
% function erb = freq2erb(freq)
%
% 1. Description:
%       Converts freq (Hz) to ERBs
% See also freqtoaud
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

erb = 21.4*log10(4.37*freq/1e3+1);


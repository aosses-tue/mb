function n = note2cent(freq1,freq2)
% function n = note2cent(freq1,freq2)
%
% 1. Description:
%
% 2. Stand-alone example:
%       freq1 = 830.6; % Gsh5
%       freq2 = 880; % A5
%       n = note2cent(freq1,freq2);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 19/09/2015
% Last update on: 19/09/2015 
% Last use on   : 19/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1200* (log10(freq2/freq1))/(log10(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function JND = Get_JND_freq_tones(f)
% function JND = Get_JND_freq_tones(f)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 27/08/2014
% Last update on: 27/08/2014 
% Last use on   : 27/08/2014 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    f = [125 250 500 1000 2000 4000 8000];
end

idx_low = find(f<=500);
idx_high = find(f>500);

JND_low = 3.6*ones( 1,length(idx_low) );
JND_high = 0.007*f(idx_high);

JND = [JND_low JND_high];

if nargout == 0
    
    % Acording to Figure 7.8, Fastl2007
    figure;
    semilogx(f,JND);
    xlabel('Frequency [Hz]')
    ylabel('JND in frequency, 2 \Delta f [Hz]')
    grid on;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

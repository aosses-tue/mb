function stats = Get_stats_from_audio(x,ti,tf,fs)
% function stats = Get_stats_from_audio(x,ti,tf,fs)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/10/2014
% Last update on: 16/10/2014 % Update this date manually
% Last use on   : 16/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    
    idx = 1:length(x);
    
else
    
    tmp = 1:length(x);
    idx = find(tmp >= ti*fs & tmp <= tf*fs); 
    
end

stats.mean  = mean(x(idx));
stats.std   = std(x(idx));

N = 100;
[y centres] = Probability_density_function(x(idx),N);

stats.pdf = y;
stats.pdf_centres = centres;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

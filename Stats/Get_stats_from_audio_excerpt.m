function [stats tbuf] = Get_stats_from_audio_excerpt(x,fs,bufsize,hopsize)
% function [stats tbuf] = Get_stats_from_audio_excerpt(x,fs,bufsize,hopsize)
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

ovsize  = bufsize - hopsize;

xbuf    = buffer(x,bufsize,ovsize,'nodelay');

t       = 1:size(x,1);
tbuf    = buffer(t,bufsize,ovsize,'nodelay'); 
tbuf = tbuf(1,:) / fs;

stats.mean  = mean(xbuf);
stats.std   = std(xbuf);

N = 100;
[y centres] = Probability_density_function(xbuf,N);

stats.pdf = y;
stats.pdf_centres = centres;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

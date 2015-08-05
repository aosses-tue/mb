function p = setupStandardParams(p)
% function p = setupStandardParams(p)
%
% Some default values for ACE and F0mod maps used by Matthias on his thesis
%
% Some parameters assign:
%       DEBUG = 0           NO debug
%       channel_stim_rate   always 1800 pps
%       analysis_rate       always 1800 pps
%
% Programmed by Matthias Milczynski
% Comments by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
   p = []; 
end
p = Ensure_field(p, 'DEBUG', 0);
p = Ensure_field(p, 'channel_stim_rate', 1.8e3);
p = Ensure_field(p, 'analysis_rate', 1.8e3);
p = Ensure_field(p, 'num_selected', 8);
p = Ensure_field(p, 'micfront', 0);
p = Ensure_field(p, 'processPreFiltered', 0);
p = Ensure_field(p, 'normalize', 1);
p = Ensure_field(p, 'filterbank', 'Freedom');
p = Ensure_field(p, 'qic', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
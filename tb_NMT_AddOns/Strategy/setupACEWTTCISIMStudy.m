function ace = setupACEWTTCISIMStudy(ace)
% function ace = setupACEWTTCISIMStudy(ace)
%
% Creates structure 'ace' with all the parameters needed for a default
% CISimulation (Vocoder)
%
% Dependencies:
%       ACEResynth_map
%
% % Example:
%       ace = setupACEWTTCISIMStudy;
%
% Programmed by Matthias Milczynski, comments by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    ace = [];
end
ace = Ensure_field(ace, 'resynth_order', 'postFilt');
ace = Ensure_field(ace, 'normalize', 0);
ace = Ensure_field(ace, 'analysis_rate', 1800);
ace = Ensure_field(ace, 'num_bands', 22);
ace = Ensure_field(ace, 'num_selected', 8);
ace = ACEResynth_map(ace); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
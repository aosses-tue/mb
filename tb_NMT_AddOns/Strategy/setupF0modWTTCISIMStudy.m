function f0mod = setupF0modWTTCISIMStudy(f0mod)
% function f0mod = setupF0modWTTCISIMStudy(f0mod)
%
% Creates structure 'f0mod' with all the parameters needed for a default
% CISimulation (Vocoder)
%
% Dependencies:
%   F0mod128Resynth_map
%
% % Example:
%       f0m = setupF0modWTTCISIMStudy;
%
% Programmed by Matthias Milczynski, comments by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files   = { 'Matlab_version.m', ...
                    'LPF_zero_ph_proc.m'};
if isunix 
    UserName = 'alejandro';
else
    UserName = 'Administrator';
end
    
Check_dependencies_ExpORL(List_of_files, UserName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    f0mod = []; 
end
f0mod = Ensure_field(f0mod, 'resynth_order', 'postFilt');
f0mod = Ensure_field(f0mod, 'normalize', 0);
f0mod = Ensure_field(f0mod, 'analysis_rate', 1800);
f0mod = Ensure_field(f0mod, 'num_bands', 22);
f0mod = Ensure_field(f0mod, 'num_selected',  8);
f0mod = Ensure_field(f0mod, 'cisim ', 1);
f0mod = Ensure_field(f0mod, 'f0File', 1);
f0mod = Ensure_field(f0mod, 'rmsequalize', 0);
f0mod = Ensure_field(f0mod, 'removeF0Peaks', 0);
f0mod = F0mod128Resynth_map(f0mod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
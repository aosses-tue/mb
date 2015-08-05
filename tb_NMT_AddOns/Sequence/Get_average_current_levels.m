function v = Get_average_current_levels(seq, idx_min, idx_max, NumMaxima)
% function v = Get_average_current_levels(seq, idx_min, idx_max, NumMaxima)
%
% Same than Get_mean_current_levels, but considers null current values.
% Returns a vector of mean current level on each electrode.
%
% v = Get_average_current_levels(seq)
%
% seq: Pulse sequence.
% v:   Column vector of current levels.
% 
% Edited by Alejandro Osses, ExpORL, 2013
% Adapted from Get_mean_current_levels, NMT by Brett Swanson (Change 86418, Cochlear Ltd.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    return;
end

if nargin < 2
    idx_min     = 1;
end

if nargin < 3
    idx_max     = length(seq.electrodes);
end

if nargin < 4
    NumMaxima   = 8;
end

Total_num_pulses = ( idx_max-idx_min )/NumMaxima;

seq.electrodes  = seq.electrodes(idx_min:idx_max);
seq.current_levels = seq.current_levels(idx_min:idx_max);

v = zeros(22, 1);
for e = 1:22
	x = (seq.electrodes == e);
    
    current_levels = seq.current_levels(x);
	if ~isempty(current_levels)
        v(e) = sum(current_levels)/Total_num_pulses;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

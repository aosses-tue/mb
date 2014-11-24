function [seq, collision_channels] = Remove_collisions(seq, min_period)

% Remove_collisions: Remove pulses that overlap in time from a chan-mag sequence.
% To resolve a collision, the pulse with smaller amplitude is deleted from the sequence.
% The sequence field periods defines the time from the start of the pulse
% to the start of the next pulse.
% A collision is defined to occur if a period is smaller than the minimum period.
% The minimum period is equal to the RF frame width plus the minimum RF frame gap.
% It is assumed all pulses have the same min_period (i.e. same phase width, phase gap).
%
% [seq_out, collision_channels] = Remove_collisions(seq_in, min_period)
%
% seq_in:             Input sequence, with pulses that may overlap in time.
% min_period:         Minimum period allowed between pulses.
% seq_out:            Output sequence, with strictly sequential pulses.
% collision_channels: Channels that had pulses removed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_stims = Get_num_pulses(seq);
[pulse_times, duration] = Get_pulse_times(seq);

k = 1;
m = 1;
collision_channels = [];
for s = 1:(num_stims-1);

	period = pulse_times(k+1) - pulse_times(k);
		
	if (period >= min_period)
		% No collision
		k = k + 1;
	else
		% Collision:
		% Delete the smaller magnitude stimulus:
		if (seq.magnitudes(k+1) < seq.magnitudes(k))
			collision_index = k+1;
		else
			collision_index = k;
		end
		collision_channels(m) = seq.channels(collision_index);
		m = m + 1;
		seq.magnitudes(collision_index) = [];
		seq.channels  (collision_index) = [];
		pulse_times   (collision_index) = [];

		% Note: k is not incremented,
		% because we have to check for further collisions.
	end
end

if (not(isempty(collision_channels)))
	collision_channels = unique(collision_channels);
end	

periods			= diff(pulse_times);
last_period		= periods(end);
seq.periods		= [periods; last_period];
% seq.periods = diff([pulse_times; duration]);

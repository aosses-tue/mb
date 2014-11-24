function [t, duration] = Get_pulse_times(seq)

% Get_pulse_times: Get the time of each pulse of a sequence.
% The sequence must have a field "periods".
% The period is the time from the start of the pulse,
% to the start of the next pulse.
% i.e. the pulse starts at the start of the period.
% The first pulse is defined to occur at time 0.
% The last time (the sequence duration) can be returned as a second output;
% note that there is no pulse there.
% This can be convenient, because:
% seq.periods == diff([t; duration])
%
% [t, duration] = Get_pulse_times(seq)
%
% seq:      Sequence struct
% seq.periods: time from start of pulse to start of next pulse.
%
% t:        Time of the start of each pulse.
% duration: Duration of the sequence (includes last period).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_pulses = Get_num_pulses(seq);

if length(seq.periods) == 1
	t = seq.periods * (0:num_pulses)';
else
	t = [0; cumsum(seq.periods)];
end

duration = t(end);
t(end) = [];
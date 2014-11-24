function q = Overlap_sequences(q1, q2, time_offset)

% Overlap_sequences: Overlap two sequences in time.
% The sequences must both have a field "periods".
% Collisions can be produced. 
%
% q = Overlap_sequences(q1, q2, time_offset)
%
% q1:          Sequence 1.
% q2:	       Sequence 2.
% time_offset: Delay from start of sequence 1 to start of sequence 2.
%
% q:           Overlapped sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (time_offset < 0)
	q = Overlap_sequences(q2, q1, -time_offset);
	return;
end

[times1, duration1] = Get_pulse_times(q1);
[times2, duration2] = Get_pulse_times(q2);

duration = max([duration1; duration2 + time_offset]);

[sorted_times, indices] = sort([times1; times2 + time_offset]);

q1r = rmfield(q1, 'periods');
q2r = rmfield(q2, 'periods');
q = Concatenate_sequences(q1r, q2r);

q.channels   = q.channels(indices);
q.magnitudes = q.magnitudes(indices);
q.periods    = diff([sorted_times; duration]);	
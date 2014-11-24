function t = Get_sequence_duration(seq)

% Get_sequence_duration: Get the duration of a sequence.
% The sequence must have a field "periods".
%
% t = Get_sequence_duration(seq)
%
% seq: Sequence struct.
% t:   Duration (in microseconds).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_pulses = Get_num_pulses(seq);

if length(seq.periods) == 1
	t = seq.periods * num_pulses;
else
	t = sum(seq.periods);
end

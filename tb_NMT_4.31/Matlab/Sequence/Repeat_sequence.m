function qo = Repeat_sequence(qi, num_copies, period)

% Repeat_sequence: Repeat a sequence periodically.
% The period is allowed to be less than the sequence duration,
% in which case the sequences will overlap,
% and collisions can occur.
%
% qo = Repeat_sequence(qi, num_copies, period)
%
% qi:         Input sequence.
% num_copies: Number of repetitions of input sequence.
% period:     Repetition period.
%
% qo:         Output sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (num_copies == 1)
	qo = qi;
	return;
end

[pulse_times, duration] = Get_pulse_times(qi);

start_times	= (0:num_copies-1) * period;
start_times	= repmat(start_times, size(pulse_times));

pulse_times	= repmat(pulse_times, 1, num_copies);
pulse_times	= pulse_times + start_times;
[sorted_times, indices] = sort(pulse_times(:));

% Copy all fields except periods:
qu = rmfield(qi, 'periods');
qr = Concatenate_sequences(qu, num_copies);

% Sort the remaining fields in time order:
qo.channels   = qr.channels(indices);
qo.magnitudes = qr.magnitudes(indices);

duration	= duration + (num_copies-1) * period;
qo.periods	= diff([sorted_times; duration]);

function [intervals, counts] = Interval_histogram(seq)

% Interval_histogram: Calculates the intervals between stim pulses in a sequence.
%
% [intervals, count] = Interval_histogram(seq)
%
% seq:	Pulse sequence.
%
% intervals: A vector of the intervals between stimulation pulses, in microsecs,
%              sorted from smallest to largest.
% counts:    A vector of the number of occurrences of each interval.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seq = Remove_idles_from_sequence(seq);

% Convert microseconds to RF cycles:
% (avoids floating point comparison issues)
num_cycles = round(5 * seq.periods);

unique_num_cycles = unique(num_cycles);

counts = zeros(size(unique_num_cycles));
for n = 1:length(unique_num_cycles)
	t = unique_num_cycles(n);
	x = (num_cycles == t);
	counts(n) = sum(x);
end

intervals = unique_num_cycles/5;	% convert back to microseconds.
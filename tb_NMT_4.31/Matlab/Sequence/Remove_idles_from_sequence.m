function seq = Remove_idles_from_sequence(seq)

% Remove_idles_from_sequence: Remove idle pulses from a sequence.
% The timing of the remaining stimulus pulses is maintained;
% a constant-period sequence will be converted to a variable-period sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seq = Complete_sequence(seq);
seq = Remove_pulses_from_sequence(seq, seq.magnitudes < 0);
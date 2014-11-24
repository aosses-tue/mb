function seq = Time_reverse_sequence(seq)

% Time_reverse_sequence: Reverses a sequence in time
%
% qo = Time_reverse_sequence(qi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seq.channels	= flipud(seq.channels);
seq.magnitudes	= flipud(seq.magnitudes);

periods_ex		= seq.periods(1:end-1);
seq.periods		= [ flipud(periods_ex); seq.periods(end)];

function seq = Trim_sequence(seq)
% Trim_sequence: Trim idle pulses from the start and end of a sequence.
% This is useful for sequences captured from the PCI by RxFrames or similar.
% function seq = Trim_sequence(seq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(1,1,nargin));

% If argument is a string, assume it is a file name:

if ischar(seq)
	seq = Read_sequence(seq);
end

stim_indices = find(seq.current_levels > 0);
seq = Subsequence(seq, stim_indices(1):stim_indices(end));


function qq = Stack_sequences(cell_seqs)

% Stack_sequences: Merge single channel sequences into one sequence.
% Each input sequence is put onto a separate channel of the output sequence.
% The main use is to allow several single channel sequences to be 
% plotted together by Plot_sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:length(cell_seqs)
	cell_seqs{n}.channels = n;
	cell_seqs{n} = Remove_idles_from_sequence(cell_seqs{n});
end

qq = cell_seqs{1};
for n = 2:length(cell_seqs)
	qq = Overlap_sequences(qq, cell_seqs{n}, 0);
end

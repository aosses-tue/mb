function seq = Subsequence(seq, indices)
% Subsequence: returns a sequence containing a subset of the pulses.
% function seq = Subsequence(seq, indices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names  = fieldnames(seq);
num_fields   = length(field_names);

for f = 1:num_fields
	name = field_names{f};
	vec  = getfield(seq, name);
	if (length(vec) > 1)
		seq = setfield(seq, name, vec(indices));
	end
end

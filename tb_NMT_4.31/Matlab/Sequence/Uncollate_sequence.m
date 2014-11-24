function v = Uncollate_sequence(q, num_bands)

% Uncollate_sequence: Convert a (channel,magnitude) sequence into a FTM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
	num_bands  = max(q.channels);
end
num_pulses = Get_num_pulses(q);
v = zeros(num_bands, num_pulses);
for n = 1:num_pulses
	if q.channels(n) > 0
		v(q.channels(n), n) = q.magnitudes(n);
	end
end

function result = Trim_sequence_test

% Trim_sequence_test: Test of Trim_sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODE = Implant_modes;

rand('state', 0);

n = [3 5 2 7 2]; 
num_pulses	= sum(n);

seq.electrodes		= ceil(22*rand(num_pulses,1));
seq.modes			= MODE.MP1;
seq.current_levels	= [	      zeros(n(1), 1);	% trim these
						repmat(100, n(2), 1);
						      zeros(n(3), 1);	% don't trim these
						repmat(150, n(4), 1);
						      zeros(n(5), 1)];	% trim these
trim = Trim_sequence(seq);

if verbose
	Disp_sequence(seq);
	Disp_sequence(trim);
end

x1 = n(1)+1;
x2 = num_pulses - n(5);

Tester(trim.electrodes,		seq.electrodes(x1:x2));
Tester(trim.modes,			seq.modes);
Tester(trim.current_levels,	seq.current_levels(x1:x2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result.

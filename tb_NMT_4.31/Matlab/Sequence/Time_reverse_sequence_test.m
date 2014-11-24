function result = Time_reverse_sequence_test

% Time_reverse_sequence_test: Test of Time_reverse_sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform period:

u.channels		= [  1;  2;  3;  4];
u.magnitudes	= [ 10; 20; 30; 40];
u.periods		= 100;

ru = Time_reverse_sequence(u);

Tester(Get_num_pulses(u),        Get_num_pulses(ru));
Tester(Get_sequence_duration(u), Get_sequence_duration(ru));

ru_.channels	= [  4;  3;  2;  1];
ru_.magnitudes	= [ 40; 30; 20; 10];
ru_.periods		= 100;

Tester(ru,ru_);

if verbose
	Plot_sequence({u,ru},{'u','ru'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable period:

v.channels		= [  1;  2;  3;  4];
v.magnitudes	= [ 10; 20; 30; 40];
v.periods		= [100;200;300;400];

rv = Time_reverse_sequence(v);

Tester(Get_num_pulses(v),        Get_num_pulses(rv));
Tester(Get_sequence_duration(v), Get_sequence_duration(rv));

rv_.channels	= [  4;  3;  2;  1];
rv_.magnitudes	= [ 40; 30; 20; 10];
rv_.periods		= [300;200;100;400];

Tester(rv,rv_);

if verbose
	Plot_sequence({v,rv},{'v','rv'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

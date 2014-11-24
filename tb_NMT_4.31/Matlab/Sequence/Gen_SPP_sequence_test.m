function result = Gen_SPP_sequence_test

% Gen_SPP_sequence_test: Test of Gen_SPP_sequence.
% Test cases for ratio carrier_rate/modulation_freq:
% - even integer
% - odd integer
% - non-integer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gen_SPP_sequence(channel, magnitude, carrier_rate, modulation_rate, duration)

% modulation_freq == carrier_rate/10
u100 = Gen_SPP_sequence(2, 9, 1000, 100, 0.08);
u100_.channels		= 2;
one_period			= [9; repmat(-1,9,1)];
u100_.magnitudes	= repmat(one_period, 8, 1);
u100_.periods		= 1000;
Tester(u100,u100_);

u125 = Gen_SPP_sequence(2, 9, 1000, 125, 0.08);
u125_.channels		= 2;
one_period			= [9; repmat(-1,7,1)];
u125_.magnitudes	= [repmat(one_period, 10, 1)];
u125_.periods		= 1000;
Tester(u125,u125_);

% non-integer ratio:
% period_pulses = 1000/150 = 6.7
u150 = Gen_SPP_sequence(2, 9, 1000, 150, 0.08);
u150_.channels		= 2;
p7					= [9; repmat(-1,6,1)];
p6					= [9; repmat(-1,5,1)];
u150_.magnitudes	= repmat([p7;p7;p6], 4, 1);
u150_.periods		= 1000;
Tester(u150,u150_);

% modulation_freq == carrier_rate/5 (odd integer)
u200 = Gen_SPP_sequence(2, 9, 1000, 200, 0.08);
u200_.channels		= 2;
one_period			= [9;-1;-1;-1;-1];
u200_.magnitudes	= repmat(one_period, 16, 1);
u200_.periods		= 1000;
Tester(u200,u200_);

% non-integer ratio:
u400 = Gen_SPP_sequence(2, 9, 1000, 400, 0.08);
u400_.channels		= 2;
one_period			= [9;-1;-1;9;-1];
u400_.magnitudes	= repmat(one_period, 16, 1);
u400_.periods		= 1000;
Tester(u400,u400_);

% modulation_freq == carrier_rate/2
u500 = Gen_SPP_sequence(2, 9, 1000, 500, 0.08);
u500_.channels		= 2;
u500_.magnitudes	= repmat([9;-1], 40, 1);
u500_.periods		= 1000;
Tester(u500,u500_);

% single pulse:
u1 = Gen_SPP_sequence(3, 30, 2500, 100, 0);
u1_.channels	= 3;
u1_.magnitudes	= 30;
u1_.periods		= 400;
Tester(u1,u1_);

if verbose > 1
	Plot_sequence({u100,u125,u150,u200,u400,u500},{'100','125','150','200','400','500'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

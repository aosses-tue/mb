function result = Gen_sequence_test

% Gen_sequence_test: Test of Gen_sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gen_sequence(channel, magnitude, stim_rate, duration)

% single pulse:
u1 = Gen_sequence(3, 30, 250, 4e-3);
u1_.channels	= 3;
u1_.magnitudes	= 30;
u1_.periods		= 4000;
Tester(u1,u1_);

% We always get at least one pulse:
u0 = Gen_sequence(3, 30, 250, 0);
Tester(u0,u1_);

% 2 pulses:
u2 = Gen_sequence(2, 0.2, 1000, 2e-3);
u2_.channels	= 2;
u2_.magnitudes	= [0.2; 0.2];
u2_.periods		= 1000;
Tester(u2,u2_);

% many pulses:
u3 = Gen_sequence(10, 100, 2000, 1);
u3_.channels	= 10;
u3_.magnitudes	= repmat(100, 2000, 1);
u3_.periods		= 500;
Tester(u3,u3_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

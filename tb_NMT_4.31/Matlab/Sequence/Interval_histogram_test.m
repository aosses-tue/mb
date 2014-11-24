function result = Interval_histogram_test

% Interval_histogram_test: Test of Interval_histogram.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform rate sequences:

qu1 = Gen_sequence(7, 1, 2000, 0.01);
[tu1,cu1] = Interval_histogram(qu1);

Tester(tu1,  500);
Tester(cu1,   20);

qu2 = Gen_sequence(6, 0.5, 1000, 0.005);
[tu2,cu2] = Interval_histogram(qu2);

Tester(tu2, 1000);
Tester(cu2,    5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate two uniform rate sequences:

q12 = Concatenate_sequences(qu1, qu2);
[t12,c12] = Interval_histogram(q12);

Tester(t12, [ 500;1000]);
Tester(c12, [  20;   5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable-rate sequence:

qv1.channels	= 5;
qv1.magnitudes	= 1;
qv1.periods		= [
	100;
	200;
	100;
	300;
	300;
	100;
	200;
	500;
	100;
	300;
	100;
	200
	];

[tv1,cv1] = Interval_histogram(qv1);

Tester(tv1, [ 100; 200; 300; 500]);
Tester(cv1, [   5;   3;   3;   1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable-rate sequence containing idle pulses:

qv2 = qv1;
% Convert some pulses to idles:
qv2.magnitudes	= [
	 1;	%		  100
	 1;	%		  200
	 1;	% 100 } = 700
	-1;	% 300 }
	-1;	% 300 }
	 1;	%		  100
	 1;	% 200 } = 700
	-1;	% 500 }
	 1;	% 		  100
	 1;	%		  300
	 1;	%		  100
	 1	%		  200
	];

[tv2,cv2] = Interval_histogram(qv2);

Tester(tv2, [ 100; 200; 300; 700]);
Tester(cv2, [   4;   2;   1;   2]);

if verbose
	Plot_sequence({qv1,qv2});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

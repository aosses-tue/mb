function result = Repeat_sequence_test

% Repeat_sequence_test: Test of Repeat_sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q.channels		= [  1;  2;  3];
q.periods		= [100;200;300];
q.magnitudes	= 100 + 10 * q.channels;

% one copy: no action:
a1 = Repeat_sequence(q, 1, 100);
Tester(a1, q);

% period > duration
a2 = Repeat_sequence(q, 2, 700);
a2_.channels	= [  1;  2;  3;  1;  2;  3];
a2_.periods		= [100;200;400;100;200;300];
a2_.magnitudes	= 100 + 10 * a2_.channels;
Tester(a2, a2_);

% period == duration
b2 = Repeat_sequence(q, 2, Get_sequence_duration(q));
b2_.channels	= [  1;  2;  3;  1;  2;  3];
b2_.periods		= [100;200;300;100;200;300];
b2_.magnitudes	= 100 + 10 * b2_.channels;
Tester(b2, b2_);

% period < duration; no overlap
c2 = Repeat_sequence(q, 2, 400);
c2_.channels	= [  1;  2;  3;  1;  2;  3];
c2_.periods		= [100;200;100;100;200;300];
c2_.magnitudes	= 100 + 10 * c2_.channels;
Tester(c2, c2_);

% period < duration; collision with last pulse
d2 = Repeat_sequence(q, 2, 300);
d2_.channels	= [  1;  2;  3;  1;  2;  3];
d2_.periods		= [100;200;  0;100;200;300];
d2_.magnitudes	= 100 + 10 * d2_.channels;
Tester(d2, d2_);

% period < duration; partial interleave
d2 = Repeat_sequence(q, 2, 250);
d2_.channels	= [  1;  2;  1;  3;  2;  3];
d2_.periods		= [100;150; 50; 50;200;300];
d2_.magnitudes	= 100 + 10 * d2_.channels;
Tester(d2, d2_);

% period < duration; complete interleave
e2 = Repeat_sequence(q, 2,  50);
e2_.channels	= [  1;  1;  2;  2;  3;  3];
e2_.periods		= [ 50; 50; 50;150; 50;300];
e2_.magnitudes	= 100 + 10 * e2_.channels;
Tester(e2, e2_);

% period < duration; partial interleave
d4 = Repeat_sequence(q, 4, 250);
d4_.channels	= [  1;  2;  1;  3;  2;  1;  3;  2;  1;  3;  2;  3];
d4_.periods		= [100;150; 50; 50;150; 50; 50;150; 50; 50;200;300];
d4_.magnitudes	= 100 + 10 * d4_.channels;
Tester(d4, d4_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose
	Plot_sequence({q,a1,a2,b2,c2,d2,e2,d4},{'q','a1','a2','b2','c2','d2','e2','d4'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

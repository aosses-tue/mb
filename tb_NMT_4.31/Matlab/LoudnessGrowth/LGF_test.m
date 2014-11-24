function result = LGF_test

% LGF_test: Test of LGF_proc, LGF_sequence_proc, LGF_inv_proc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax = 169;
x = 0:xmax;

Qs = 10:10:50;
Bs =  0: 4:16;

v = zeros(length(x), length(Bs), length(Qs));

for b = 1:length(Bs)
	for q = 1:length(Qs)
	
		p = [];
		p.base_level = Bs(b);
		p.Q          = Qs(q);
		p.sat_level  = 150;
		p.sub_mag    = 0;		

		p = LGF_proc(p);
		
		y = LGF_proc(p, x);
		
		v(:,b,q) = y;	% save in a 3D array

		% Make some simple checks on the output:
		
		% All inputs <= base_level should give output = 0:
		sub_inputs  = (x <= p.base_level);
		sub_outputs = y(sub_inputs);
		Tester(sub_outputs == 0);
		
		% All inputs > sat_level should give output = 1:
		sat_inputs  = (x >= p.sat_level);
		sat_outputs = y(sat_inputs);
		Tester(sat_outputs == 1);
		
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expected = load('LGF_test.mat');	% read manually-checked result.
Tester(v, expected.v, 1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse processing:

p.base_level =   4/256;
p.sat_level  = 150/256;
p = LGF_proc(p);

u = (0:256)/256;
v = LGF_proc(p, u);

p = LGF_inv_proc(p);
w = LGF_inv_proc(p, v);

expected_equal_range = (u >= p.base_level) & (u <= p.sat_level);
Tester(u(expected_equal_range), w(expected_equal_range), 1e-10);

if verbose >= 2
	figure
	hold on
	plot(u, v, 'g');
	plot(w, v, 'r');

	plot(u, w, 'c');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence processing:

qi.channels = 2;
qi.magnitudes = u;
qo = LGF_sequence_proc(p, qi);
Tester(qo.channels, qi.channels);
Tester(qo.magnitudes, v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

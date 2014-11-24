function result = LGF_alpha_test

% LGF_alpha_test: Test of LGF_alpha
% Exercise LGF_alpha for all standard Base-level and Q values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Colin Irwin, Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qs = 10:50;
Bs =  1:15;
sat_level = 150;

for b = 1:length(Bs);
	B = Bs(b);
	for q = 1:length(Qs);
		Q = Qs(q);
		alpha(b,q) = LGF_alpha(Q, B, sat_level);
		Qdiff(b,q)  = Q - LGF_Q(alpha(b,q), B, sat_level);
	end
end

max_Qdiff = max(max(abs(Qdiff)));

if verbose;
	disp(sprintf('Worst error = %g', max_Qdiff));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tester(max_Qdiff < 1e-6);
result = Tester;	% Report result

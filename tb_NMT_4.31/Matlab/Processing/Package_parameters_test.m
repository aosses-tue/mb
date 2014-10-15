function result = Package_parameters_test

% Package_parameters_test: Test Package_parameters_proc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copy1right: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normal scheme, without Package_parameters_proc:
p1  = [];
p1  = Append_process(p1, 'Example1_proc');
p1  = Append_process(p1, 'Example2_proc');
p1_.f1 = 10;
p1_.f2 = 20;
p1_.f3 = 30;
p1_.g1 = 2;
p1_.g2 = 4;
p1_.g3 = 8;
p1_.processes = {'Example1_proc'; 'Example2_proc'};
Tester(p1,  p1_);

y1 = Process(p1, 5);
Tester(y1, 158);

% Call Package_parameters_proc directly (not the usual case):
% Should return parameter struct unmodified:
p2 = Package_parameters_proc(p1);
Tester(p1,  p2);
% Processing:
% Output should be param struct, with additional field signal, equal to output:
y2 = Package_parameters_proc(p2, 5);
y2_ = p1_;
y2_.signal = 158;

% Append Package_parameters_proc:
% Should simply append itself to processes:
p3 = Append_process(p1, 'Package_parameters_proc');
p3_ = p1_;
p3_.processes = {'Example1_proc'; 'Example2_proc'; 'Package_parameters_proc'};
Tester(p3,  p3_);
% Processing:
% Output should be param struct, with additional field signal, equal to output:
y3 = Process(p3, 5);
y3_ = p3_;
y3_.signal = 158;
Tester(y3,  y3_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result.

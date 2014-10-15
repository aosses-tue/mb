function result = Split_process_test

% Split_process_test: Test Split_process.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Cochlear Ltd
%   $Header: //dsp/Nucleus/_Releases/NMT_4.30/Matlab/Processing/Split_process_test.m#1 $
%   $Change: 86418 $
% $DateTime: 2008/03/04 14:27:13 $
%   Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qc  = [];
qc  = Append_process(qc, 'Example1_proc');
qc  = Append_process(qc, 'Example2_proc');
qc  = Append_process(qc, 'Abs_proc');

% Remove processes from head:
[qc_head0, qc_tail0] = Split_process(qc, 0);

qc_head0_ = qc;
qc_head0_.processes = {};

Tester(qc_head0, qc_head0_);
Tester(qc_tail0, qc);

% Split after 1st process:
[qc_head1, qc_tail1] = Split_process(qc, 1);

qc_head1_ = qc;
qc_tail1_ = qc;
qc_head1_.processes = {'Example1_proc'};
qc_tail1_.processes = {'Example2_proc'; 'Abs_proc'};

Tester(qc_head1, qc_head1_);
Tester(qc_tail1, qc_tail1_);

% Split after 2nd process:
[qc_head2, qc_tail2] = Split_process(qc, 2);

qc_head2_ = qc;
qc_tail2_ = qc;
qc_head2_.processes = {'Example1_proc'; 'Example2_proc'};
qc_tail2_.processes = {'Abs_proc'};

Tester(qc_head2, qc_head2_);
Tester(qc_tail2, qc_tail2_);

% Split after 3rd process (remove processes from tail):
[qc_head3, qc_tail3] = Split_process(qc, 3);

qc_head3_ = qc;
qc_tail3_ = qc;
qc_tail3_.processes = {};

Tester(qc_head3, qc_head3_);
Tester(qc_tail3, qc_tail3_);

% Split before last process:
% (same as after 2nd process)
[qc_head_m1, qc_tail_m1] = Split_process(qc, -1);

Tester(qc_head_m1, qc_head2_);
Tester(qc_tail_m1, qc_tail2_);

% Split before 2nd last process:
% (same as after 1st process)
[qc_head_m2, qc_tail_m2] = Split_process(qc, -2);

Tester(qc_head_m2, qc_head1_);
Tester(qc_tail_m2, qc_tail1_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split process by name:

% Split after 1st process, by name:
[qc_head1, qc_tail1] = Split_process(qc, 'Example1_proc');
Tester(qc_head1, qc_head1_);
Tester(qc_tail1, qc_tail1_);

% Split after 1st process, by partial name:
[qc_head1, qc_tail1] = Split_process(qc, 'Example1');
Tester(qc_head1, qc_head1_);
Tester(qc_tail1, qc_tail1_);

% Split after 2nd process, by name:
[qc_head2, qc_tail2] = Split_process(qc, 'Example2_proc');
Tester(qc_head2, qc_head2_);
Tester(qc_tail2, qc_tail2_);

% Split after 2nd process, by partial name:
[qc_head2, qc_tail2] = Split_process(qc, 'Example2');
Tester(qc_head2, qc_head2_);
Tester(qc_tail2, qc_tail2_);

% Partial name finds first match:
[qc_head1, qc_tail1] = Split_process(qc, 'Example');
Tester(qc_head1, qc_head1_);
Tester(qc_tail1, qc_tail1_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result.

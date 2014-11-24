function result = Blank_sequence_test;

% Blank_sequence_test: Test of Blank_sequence_proc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Tim Neal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1000;
p.implant_stim_rate = 1;
x.electrodes = ones(N,1);
x.current_levels = ones(N,1);
x.periods = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1: Blank-out entire sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.blank_duration = N;
p.blank_freq = 1/N;

p = Blank_sequence_proc(p);
y = Blank_sequence_proc(p, x);

output = zeros(N,1);
Tester(y.electrodes,output);
Tester(y.current_levels,output);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2: Blank-out half the sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.blank_duration = N/2;
p.blank_freq = 1/N;

p = Blank_sequence_proc(p);
y = Blank_sequence_proc(p, x);

output = [zeros(N/2,1); ones(N/2,1)];
Tester(y.electrodes,output);
Tester(y.current_levels,output);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3: 25% blanking duty cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.blank_duration = N/32;
p.blank_freq = 8/N;

p = Blank_sequence_proc(p);
y = Blank_sequence_proc(p, x);

output = repmat([zeros(floor(N/32)+1,1); ones(floor(3*N/32),1)],8,1);
Tester(y.electrodes,output);
Tester(y.current_levels,output);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collate results from tests...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = Tester;

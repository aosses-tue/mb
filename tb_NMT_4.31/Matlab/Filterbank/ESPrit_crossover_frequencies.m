function f = ESPrit_crossover_frequencies(table_num, upper_freq)

% ESPrit_crossover_frequencies: Nominal crossover frequencies of ESPrit filterbank.
% The ESPrit family of processors all have 22 filters,
% but they are scaled so that only 20 are useful at any one time.
%
% f = ESPrit_crossover_frequencies(table_num, upper_freq)
%
% table_num:  Table number.
% upper_freq: Upper frequency boundary (Hz). Filters above this are discarded.
% f:          Vector of crossover frequencies (Hz).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f9 = [
     40, %   0
     80, %   1
    160, %   2
    335, %   3
    541, %   4
    743, %   5
    945, %   6
   1146, %   7
   1346, %   8
   1547, %   9
   1768, %  10
   2031, %  11
   2333, %  12
   2680, %  13
   3078, %  14
   3572, %  15
   4185, %  16
   4903, %  17
   5744, %  18
   6730, %  19
   7885, %  20
   9238, %  21
  10823  %  22
];

if nargin < 1
	table_num = 9;
end
if nargin < 2
	upper_freq = f9(end);
end

% Clock division factor of ESPrit Switched Capacitor Filterbank:
n = 17 - table_num;
r = 8 ./ n;
f = f9 * r;

enabled = (f >= 100) & (f <= upper_freq);
f = f(enabled);
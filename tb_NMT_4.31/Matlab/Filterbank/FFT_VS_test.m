function result = FFT_VS_test(update_master)

% FFT_VS_test: Regression test of FFT_filterbank_proc & Vector_sum_proc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

audio = wavread('asa');

p1 = [];
p1 = Append_process(p1, 'FFT_filterbank_proc');
p1 = Append_process(p1, 'Vector_sum_proc');

v1 = Process(p1, audio);
if verbose > 2
	GUI_FTM(p1, v1, 'Filterbank output');
end

% If behaviour ever changes, we need to update master file:
if nargin > 0	
	save FFT_VS_test.mat v1
end

master = load('FFT_VS_test');
Tester(master.v1, v1, eps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

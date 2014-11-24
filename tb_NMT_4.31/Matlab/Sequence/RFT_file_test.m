function result = RFT_file_test

% RFT_file_test: Test of Read_sequence for RFT files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_mp12 = Read_sequence('mp12.rft');
Open_log('q_mp12.log', 'new');
Log_sequence(q_mp12);
expected_log = {
'     electrodes          modes current_levels   phase_widths     phase_gaps        periods'
'             11          MP1+2            110           99.8             15            499'
'             12                           120           99.8             15            499'
'             13                           130           99.8             15            499'
'             14                           140           99.8             15            499'
'              0                             0           99.8             15            499'
'             16                           160           99.8             15            499'
};
Tester(Compare_log_file(expected_log));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

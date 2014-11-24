function result = ACE_test

% ACE_test: Regression test of ACE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The gains are adjusted to match results from version 1.0.
% In version 1.0, the weights were not equalised,
% and a gain of 1/40 was applied.

asa = wavread('asa.wav');
tapestry = wavread('tapestry.wav');

p1 = [];
p1.num_bands = 22;
p1.num_selected = 12;
p1.analysis_rate = 1000;
p1.threshold_levels = repmat(  0, p1.num_bands, 1);
p1.comfort_levels   = repmat(255, p1.num_bands, 1);
p1 = ACE_map(p1);
p1.gains_dB = To_dB(sqrt(p1.power_gains) / 40);

ACE_12_1000_asa  = Process(p1, asa);
ACE_12_1000_asa_ = Read_sequence('ACE_12_1000_asa.quf');
Tester(ACE_12_1000_asa, ACE_12_1000_asa_);

p2 = [];
p2.num_bands = 6;
p2.analysis_rate = 2000;
p2.electrodes = (6:-1:1)';
p2.threshold_levels = repmat(  0, p2.num_bands, 1);
p2.comfort_levels   = repmat(255, p2.num_bands, 1);
p2 = CIS_map(p2);
p2.gains_dB = To_dB(sqrt(p2.power_gains) / 40);

CIS_06_2000_asa  = Process(p2, asa);
CIS_06_2000_asa_ = Read_sequence('CIS_06_2000_asa.quf');
Tester(CIS_06_2000_asa, CIS_06_2000_asa_);

p3 = [];
p3.num_bands = 22;
p3.analysis_rate = 500;
p3.threshold_levels = repmat(  0, p3.num_bands, 1);
p3.comfort_levels   = repmat(255, p3.num_bands, 1);
p3 = CIS_map(p3);
p3.gains_dB = To_dB(sqrt(p3.power_gains) / 40);

CIS_22_0500_asa  = Process(p3, asa);
CIS_22_0500_asa_ = Read_sequence('CIS_22_0500_asa.quf');
Tester(CIS_22_0500_asa, CIS_22_0500_asa_);

CIS_22_0500_tapestry  = Process(p3, tapestry);
CIS_22_0500_tapestry_ = Read_sequence('CIS_22_0500_tapestry.quf');
Tester(CIS_22_0500_tapestry, CIS_22_0500_tapestry_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

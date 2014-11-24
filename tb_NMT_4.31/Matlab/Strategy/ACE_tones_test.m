function result = ACE_tones_test

% ACE_tones_test: Test of ACE with pure tones.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pure tones:
% To measure steady-state response to pure tones,
% the length of the signal can be the same as the FFT length,
% and the analysis rate can be chosen to have no overlap between blocks,
% so the FTM will be a single column.
%
% The gains are adjusted to match results from version 1.0.
% In version 1.0, the weights were not equalised,
% and a gain of 1/40 was applied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 128;

% Common parameters:

p.audio_sample_rate	= 16000;	% Hz
p.analysis_rate		= p.audio_sample_rate/N;
p.num_bands			= 22;
p.threshold_levels	= repmat(  0, p.num_bands, 1);
p.comfort_levels	= repmat(255, p.num_bands, 1);

pt1 = p;
pt1.num_selected	=  6;
pt1 = ACE_map(pt1);
pt1.gains_dB		= To_dB(sqrt(pt1.power_gains) / 40);

n = (0:N-1)';
sin1000  = sin(2*pi*n* 1000/pt1.audio_sample_rate);			
sin2000  = sin(2*pi*n* 2000/pt1.audio_sample_rate);			
sin5000  = sin(2*pi*n* 5000/pt1.audio_sample_rate);			
sin5125  = sin(2*pi*n* 5125/pt1.audio_sample_rate);			

ACE_sin1000 = Process(pt1, sin1000);
ACE_sin2000 = Process(pt1, sin2000);
ACE_sin5000 = Process(pt1, sin5000);
ACE_sin5125 = Process(pt1, sin5125);

% Compare CIS with different idle_types:

pt2 = p;
pt2 = CIS_map(pt2);
pt2.gains_dB		= To_dB(sqrt(pt2.power_gains) / 40);

pt3 = p;
pt3.sub_mag = 0;	% values below base-level give T-level
pt3 = CIS_map(pt3);
pt3.gains_dB		= To_dB(sqrt(pt3.power_gains) / 40);

CISi_sin1000 = Process(pt2, sin1000);
CISn_sin1000 = Process(pt3, sin1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose > 1
	Disp_sequence(ACE_sin1000);
	Disp_sequence(ACE_sin2000);
	Disp_sequence(ACE_sin5000);
	Disp_sequence(ACE_sin5125);

	Disp_sequence(CISi_sin1000);
	Disp_sequence(CISn_sin1000);
end
	
% All sequences should have num_selected pulses:
Tester(Get_num_pulses(ACE_sin1000), pt1.num_selected);
Tester(Get_num_pulses(ACE_sin2000), pt1.num_selected);
Tester(Get_num_pulses(ACE_sin5000), pt1.num_selected);
Tester(Get_num_pulses(ACE_sin5125), pt1.num_selected);
	
% 5000 Hz and 5125 Hz give same envelopes, within noise floor,
% but which of the approximately equal "noise" values are rejected is unpredictable,
% so the order of the idle pulses in the sequences will differ.
pulses_sin5000 = [ACE_sin5000.electrodes, ACE_sin5000.current_levels];
pulses_sin5125 = [ACE_sin5125.electrodes, ACE_sin5125.current_levels];
Tester(sortrows(pulses_sin5000), sortrows(pulses_sin5125));

% Regression checks:
ACE_sin1000_old  = Read_sequence('ACE_sin1000.quf');
ACE_sin2000_old  = Read_sequence('ACE_sin2000.quf');
ACE_sin5000_old  = Read_sequence('ACE_sin5000.quf');
CISi_sin1000_old = Read_sequence('CISi_sin1000.quf');
CISn_sin1000_old = Read_sequence('CISn_sin1000.quf');

Tester(ACE_sin1000,  ACE_sin1000_old);
Tester(ACE_sin2000,  ACE_sin2000_old);
Tester(ACE_sin5000,  ACE_sin5000_old);
Tester(CISi_sin1000, CISi_sin1000_old);
Tester(CISn_sin1000, CISn_sin1000_old);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result


% ACE_CIS_GUI_FTM_demo: Shows ACE & CIS FTMs for a sentence.
% GUI_FTM allows audio to be reconstructed: use Audio menu.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pa = [];

% Audio sample:
file_name = 'L001s04r';
[audio, pa.audio_sample_rate] = wavread(file_name);

pa.source = 'HS8';

spk06 = pa;
spk06.num_selected	=   6;
spk06.num_bands		=  22;
spk06.analysis_rate	= 250;
spk06 = ACE_map(spk06);
spk06 = Split_process(spk06, 'Reject_smallest_proc');

ace12 = pa;
ace12.num_selected	=  12;
ace12.num_bands		=  22;
ace12.analysis_rate	= 500;
ace12 = ACE_map(ace12);
ace12 = Split_process(ace12, 'Reject_smallest_proc');

cis12 = pa;
cis12.num_bands		=   12;
cis12.analysis_rate	= 1000;
cis12 = CIS_map(cis12);
cis12 = Split_process(cis12, 'Power_sum_envelope_proc');

cis06 = pa;
cis06.num_bands		=   6;
cis06.analysis_rate	= 1000;
cis06 = CIS_map(cis06);
cis06 = Split_process(cis06, 'Power_sum_envelope_proc');

% Process with each map:

v_spk06 = Process(spk06, audio);
v_ace12 = Process(ace12, audio);
v_cis12 = Process(cis12, audio);
v_cis06 = Process(cis06, audio);

GUI_FTM(spk06, v_spk06, 'SPEAK 6/22');
GUI_FTM(ace12, v_ace12, 'ACE 12/22');
GUI_FTM(cis12, v_cis12, 'CIS 12');
GUI_FTM(cis06, v_cis06, 'CIS 6');

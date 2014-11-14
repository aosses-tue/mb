% Resample_FTM_demo: Shows effect of analysis rate lower than stim rate on ACE filterbank

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 'asa';

p.num_bands = 22;
p.channel_stim_rate = 16000/9;	% approx 1800 Hz
p.num_selected = 1;
p = Append_process(p, 'Wav_proc');
p = Append_process(p, 'FFT_filterbank_proc');
p = Append_process(p, 'Power_sum_envelope_proc');
p = Append_process(p, 'LGF_proc');
p = Append_process(p, 'Resample_FTM_proc');

p09 = p;
p09.analysis_rate = p.channel_stim_rate;
[v09, p09] = Process(p09, a);

p18 = p;
p18.analysis_rate = p.channel_stim_rate/2;
[v18, p18] = Process(p18, a);
v18(:,end) = [];

p21 = p;
p21.analysis_rate = 16000/21;
[v21, p21] = Process(p21, a);
v21(:,end) = [];

p27 = p;
p27.analysis_rate = 16000/27;
[v27, p27] = Process(p27, a);
v27(:,end-1:end) = [];

GUI_FTM(p09, {v09,v18,v21,v27}, {'1800 Hz', ' 900 Hz', ' 760 Hz', ' 600 Hz'});
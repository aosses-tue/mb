% ACE_demo: Demonstrates ACE processing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('audio') ~= 1
	audio = wavread('asa');
end

p = [];
p.num_bands = 22;
p.audio_sample_rate = 16000;
p.analysis_rate = p.audio_sample_rate;
p.num_selected  = 9;
p.channel_stim_rate = p.audio_sample_rate/10;
p.base_level    = 0;

p.source = 'HS8';

p = Append_front_end_processes(p);
p = Append_process(p, 'FFT_filterbank_proc');
p = Append_process(p, 'Vector_sum_proc');
p = Append_process(p, 'Abs_proc');
p = Append_process(p, 'LGF_proc');
p = Append_process(p, 'Uniform_sampler_proc');
p = Append_process(p, 'Reject_smallest_proc');
p = Append_process(p, 'Resample_FTM_proc');
p = Append_process(p, 'Collate_into_sequence_proc');

c = Process_chain(p, audio);
q = c{end};
Plot_sequence(q);

num_pulses = Get_num_pulses(q);
v = zeros(p.num_bands, num_pulses);
for n = 1:num_pulses
	v(q.channels(n), n) = q.magnitudes(n);
end

kk = 8:10;
GUI_FTM(p, [c(kk);{v}], [p.processes(kk); {'Sequence'}]);
%colormap(gray);


% Channel_mapping_demo: Show execution time of Channel_mapping in ACE_map.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

audio_name = 'tapestry';
audio = wavread(audio_name);

p = [];
p.num_bands = 22;
p.num_selected = 12;
p.analysis_rate = 1000;
p.threshold_levels = repmat(100, p.num_bands, 1);
p.comfort_levels   = repmat(200, p.num_bands, 1);
p = ACE_map(p);

[p, pch] = Split_process(p, -1);

tic;
chan_mag_seq = Process(p, audio);
t = toc;
fprintf('Others:          %4.2f\n', t);

tic;
pulse_seq = Process(pch, chan_mag_seq);
t = toc;
fprintf('Channel_mapping: %4.2f\n', t);

if Verbose
	Plot_sequence(pulse_seq, audio_name);
end
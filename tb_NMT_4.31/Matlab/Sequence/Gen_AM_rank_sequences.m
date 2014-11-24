function x = Gen_AM_rank_sequences(p, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = Ensure_field(x, 'rise_fraction', 0);
x.audio_sample_rate = p.channel_stim_rate;
x = Gen_pitch_rank_tones([], x);

x = Ensure_field(x, 'channel', 1);

q = Gen_sequence(x.channel, 1, p.channel_stim_rate, x.stimulus_duration);

for n = 1:x.num_stimuli
	q.magnitudes = 0.5 * (1 + x.stimuli{n});
	x.stimuli{n} = q;
end

x.silence  = Gen_sequence(1, -1, p.channel_stim_rate, x.stimulus_duration);
x.gap      = Gen_sequence(1, -1, p.channel_stim_rate, x.gap_duration);

x.stimulus_description = sprintf('AM Ch %d', x.channel);

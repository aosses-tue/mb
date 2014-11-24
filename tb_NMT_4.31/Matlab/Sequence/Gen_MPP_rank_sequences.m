function x = Gen_MPP_rank_sequences(p, x)

% Gen_MPP_rank_sequences: Generate a set of single-channel MPP sequences for pitch ranking.
% Each burst has the channel stimulation rate specified in the map.
% All sequences are on the same channel.
%
% u = Gen_MPP_rank_sequences(p, x)
%
% p:                   Map parameter struct. The only field used is:
% p.channel_stim_rate:  Channel stimulation rate.
% x:                   Experimental parameters struct.
% u:                   Output struct. A copy of x, with additional fields:
% u.stimuli:            Cell array of stimuli (channel-magnitude sequences).
% u.silence:            Silent (idle) channel-magnitude sequence, same duration as stimuli.
% u.gap:                Silent (idle) channel-magnitude sequence.
% u.stimulus_names:     A cell array of stimulus names (on,off lengths).
% u.stimulus_label:     A string to label the set of stimulus names in a GUI.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = Ensure_field(x, 'channel', 1);
x = Ensure_field(x, 'stimulus_duration',    0.50);
x = Ensure_field(x, 'gap_duration',         0.25);

x.stimuli = {};
for n = 1:size(x.burst_lengths,1)
								   %(channel, magnitude, carrier_rate, num_pulses_on, num_pulses_off, duration)
	x.stimuli       {n} = Gen_MPP_sequence(x.channel, 1, p.channel_stim_rate, x.burst_lengths(n,1), x.burst_lengths(n,2), x.stimulus_duration);
	x.stimulus_names{n} = sprintf('%2d:%2d', x.burst_lengths(n,1), x.burst_lengths(n,2));
end

x.silence = Gen_sequence(1, -1, p.channel_stim_rate, x.stimulus_duration);
x.gap     = Gen_sequence(1, -1, p.channel_stim_rate, x.gap_duration);

x.stimulus_label       = 'On:Off';
x.stimulus_description = sprintf('MPP Ch %d', x.channel);

x.period_multiples = sum(x.burst_lengths, 2);
x.freqs = p.channel_stim_rate ./ x.period_multiples;

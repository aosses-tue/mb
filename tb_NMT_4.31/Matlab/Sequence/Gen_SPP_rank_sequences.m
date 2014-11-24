function x = Gen_SPP_rank_sequences(p, x)

% Gen_SPP_rank_sequences: Generate a set of single-channel SPP sequences for pitch ranking.
% Each sequence has a different stimulation rate (the pitch ranking frequency).
% All sequences are on the same channel.
%
% u = Gen_SPP_rank_sequences(p, x)
%
% p:                     Map parameter struct. The only fields used are:
% p.implant_stim_rate:    Implant stimulation rate.
% x:                     Experimental parameters struct.
% u:                     Output struct. A copy of x, with additional fields:
% u.stimuli:              Cell array of stimuli (channel-magnitude sequences).
% u.silence:              Silent (idle) channel-magnitude sequence, same duration as stimuli.
% u.gap:                  Silent (idle) channel-magnitude sequence.
% u.stimulus_names:       A cell array of stimulus names (modulation freqs).
% u.stimulus_label:       A string to label the set of stimulus names in a GUI.
% u.stimulus_description: A string that describes the stimuli. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = Gen_pitch_rank_freqs(x);

x = Ensure_field(x, 'channel', 1);

x.stimuli = {};
for n = 1:x.num_stimuli
								  %(channel, magnitude, carrier_rate, modulation_rate, duration)
	x.stimuli{n} = Gen_SPP_sequence(x.channel, 1, p.implant_stim_rate, x.freqs(n), x.stimulus_duration);
end

x.silence  = Gen_sequence(1, -1, p.implant_stim_rate, x.stimulus_duration);
x.gap      = Gen_sequence(1, -1, p.implant_stim_rate, x.gap_duration);

x.stimulus_description = sprintf('SPP Ch %d', x.channel);

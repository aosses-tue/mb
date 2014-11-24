function x = Gen_channel_rank_sequences(p, x)

% Gen_channel_rank_sequences: Generate a set of single-channel sequences for pitch ranking.
% One channel-magnitude sequence is generated for each channel in the map.
% Each sequence has the channel stimulation rate specified in the map.
% The default behaviour is:
% - for high rates (over 400 Hz) there are no idle pulses,
% - for low rates idle pulses are inserted, so that the overall implant stimulation rate
%   is the same as specified in the map.
% This default behaviour can be overridden by the parameter x.insert_idles.
%
% u = Gen_channel_rank_sequences(p, x)
%
% p:                    Map parameter struct. The only fields used are:
% p.num_bands:            Number of bands.
% p.channel_stim_rate:    Channel stimulation rate.
% p.implant_stim_rate:    Overall implant stimulation rate.
% x:                    Experimental parameters struct.
% x.stimulus_duration:    Duration of stimulus sequences.
% x.gap_duration:         Duration of gaps between sequences.
% x.insert_idles:         Determines if idle pulses should be inserted.
% u:                    Output struct. A copy of x, with additional fields:
% u.stimuli:              Cell array of stimuli (channel-magnitude sequences).
% u.silence:              Silent (idle) channel-magnitude sequence, same duration as stimuli.
% u.gap:                  Silent (idle) channel-magnitude sequence.
% u.stimulus_names:       A cell array of stimulus names ('Ch n').
% u.stimulus_label:       A string to label the set of stimulus names in a GUI.
% u.stimulus_description: A string to describe the type of stimuli.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = Ensure_field(x, 'stimulus_duration',    0.50);
x = Ensure_field(x, 'gap_duration',         0.25);
if p.channel_stim_rate <= 400
	x = Ensure_field(x, 'insert_idles', 1);
else
	x = Ensure_field(x, 'insert_idles', 0);
end

x.stimulus_description = sprintf('%4d', round(p.channel_stim_rate));
x.stimulus_label = 'Channel';
for n = 1:p.num_bands
	x.stimulus_names{n} = sprintf('Ch %2d', n);
end

if x.insert_idles == 0
	for n = 1:p.num_bands
		x.stimuli{n} = Gen_sequence(n, 1, p.channel_stim_rate, x.stimulus_duration);
	end

	x.silence = Gen_sequence(1, -1, p.channel_stim_rate, x.stimulus_duration);
	x.gap     = Gen_sequence(1, -1, p.channel_stim_rate, x.gap_duration);
	
else
	for n = 1:p.num_bands
									  %(chan, mag, carrier_rate,        modulation_rate,     duration)
		x.stimuli{n} = Gen_SPP_sequence(n,    1,   p.implant_stim_rate, p.channel_stim_rate, x.stimulus_duration);
	end

	x.silence = Gen_sequence(1, -1, p.implant_stim_rate, x.stimulus_duration);
	x.gap     = Gen_sequence(1, -1, p.implant_stim_rate, x.gap_duration);

end


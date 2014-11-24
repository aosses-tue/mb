function seq = Gen_sequence(channel, magnitude, stim_rate, duration)

% Gen_sequence: Generate a uniform rate sequence on a single channel.
%
% seq = Gen_sequence(channel, magnitude, stim_rate, duration)
%
% channel:   Channel number for stimulation.
% magnitude: Magnitude of stimulation (constant).
% stim_rate: Channel stimulation rate (in Hertz).
% duration:  Duration of sequence (in seconds).
%
% seq:       Channel-magnitude sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantise stim_rate to RF frequency
period = round(5e6 / stim_rate) / 5;		% microseconds

num_pulses = round(1e6 * duration / period);
if (num_pulses < 1)
	num_pulses = 1;
end

seq.channels	= channel;
seq.magnitudes	= repmat(magnitude, num_pulses, 1);
seq.periods		= period;	

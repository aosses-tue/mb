function seq = Gen_MPP_sequence(channel, magnitude, carrier_rate, num_pulses_on, num_pulses_off, duration)

% Gen_MPP_sequence: Generate a Multiple-Pulse-per-Period sequence on a single channel.
%
% seq = Gen_MPP_sequence(channel, magnitude, carrier_rate, num_pulses_on, num_pulses_off, duration)
%
% channel:        Channel number for stimulation.
% magnitude:      Magnitude of stimulation (constant).
% carrier_rate:   Pulse rate; channel stimulation rate during "on" burst (in Hertz).
% num_pulses_on:  Number of pulses in each "on" burst.
% num_pulses_off: Number of pulses in each "off" burst.
% duration:       Duration of sequence (in seconds).
%
% seq:            Channel-magnitude sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create fields in standard order (channels,magnitudes,periods):
seq.channels	= channel;

% Create one period ("on" burst and "off" burst):
seq.magnitudes	= [ repmat(magnitude, num_pulses_on,  1);
					repmat(-1,        num_pulses_off, 1)];

% Quantise stim rate to RF frequency:
seq.periods		= round(5e6 / carrier_rate) / 5;		% microseconds;	

% Repeat this sequence:
modulation_period = seq.periods * (num_pulses_on + num_pulses_off);	% microseconds
num_periods = round(duration * 1e6 / modulation_period); 
seq = Concatenate_sequences(seq, num_periods);
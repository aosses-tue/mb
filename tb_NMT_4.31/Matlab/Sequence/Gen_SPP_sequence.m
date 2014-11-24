function seq = Gen_SPP_sequence(channel, magnitude, carrier_rate, modulation_rate, duration)

% Gen_SPP_sequence: Generate a single pulse per period sequence on a single channel.
%
% seq = Gen_SPP_sequence(channel, magnitude, carrier_rate, modulation_rate, duration)
%
% channel:         Channel number for stimulation.
% magnitude:       Magnitude of stimulation pulses.
% carrier_rate:    Total stimulation rate (in Hertz), including idle pulses.
% modulation_rate: Modulation rate (in Hertz).
% duration:        Duration of sequence (in seconds).
%
% seq:             Channel-magnitude sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a constant rate idle sequence:
seq = Gen_sequence(channel, -1, carrier_rate, duration);

% Rates are quantised so that all periods are integral number of RF cycles:
carrier_period_cycles    = round(5 * seq.periods);
modulation_period_cycles = round(5e6/modulation_rate);

% Change some idle pulses to stimulus pulses:
num_pulses = length(seq.magnitudes);
time_until_pulse = 0;
for n = 1:num_pulses
	if time_until_pulse <= 0
		seq.magnitudes(n) = magnitude;
		time_until_pulse  = time_until_pulse + modulation_period_cycles;
	end
	time_until_pulse = time_until_pulse - carrier_period_cycles;
end	

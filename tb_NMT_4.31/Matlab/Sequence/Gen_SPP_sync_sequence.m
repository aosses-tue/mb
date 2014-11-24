function seq = Gen_SPP_sync_sequence(channel, magnitude, carrier_rate, modulation_freq, duration)

% Gen_SPP_sync_sequence: Generate a single pulse per period sequence on a single channel.
% The sequence will be "synchronised", i.e. exactly periodic.
%
% seq = Gen_SPP_sync_sequence(channel, magnitude, carrier_rate, modulation_freq, duration)
%
% channel:         Channel number for stimulation.
% magnitude:       Magnitude of stimulation pulses.
% carrier_rate:    Total stimulation rate (in Hertz), including idle pulses.
% modulation_freq: Modulation frequency (in Hertz).
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

% Carrier rate is quantised to an integral number of RF cycles:
carrier_rate = 1e6 / seq.periods;
% Quantise modulation frequency to an integer sub-multiple of carrier rate:
modulation_period_pulses = round(carrier_rate/modulation_freq);

% Change some idle pulses to stimulus pulses:
seq.magnitudes(1:modulation_period_pulses:end) = magnitude;
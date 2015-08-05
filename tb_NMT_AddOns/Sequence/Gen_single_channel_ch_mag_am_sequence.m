function ch_mag = Gen_single_channel_ch_mag_am_sequence(duration, electrode, current_level, channel_stim_rate, mode, phase_width, phase_gap, num_selected, fm)
% function ch_mag = Gen_single_channel_ch_mag_am_sequence(duration, electrode, current_level, channel_stim_rate, mode, phase_width, phase_gap, num_selected, fm)
%
% Parameters are:
%   duration: stimulus duration in seconds
%   channel: electrode number 
%   current_level: current level in % dyn. range
%   channel_stim_rate: number of pulses per second per channel
%   num_selected: number of active channels
%   phase_gap: gap between pulses within biphasic pulse in microseconds
%   phase_width: pulse duration of single pulse in microseconds
%
% Comments modified by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ch_mag             = [];
num_pulses      = round(duration*channel_stim_rate);
single_pulse    = zeros(num_selected, 1);
period          = 1e6/(num_selected*channel_stim_rate);

if( electrode ~= 0 )
    ch = 23 - electrode;
else
    ch = 0;
end
t = 0:1/channel_stim_rate:duration - 1/channel_stim_rate;
if fm == 0
    sig = ones(length(t),1);
else
    sig = sin(2*pi*t*fm);
end
modsig = (0.5*(sig + 1))*current_level;
ch_seq = [ ch; single_pulse(2:end)]; % swap channel and electrode

cl_seq = zeros(num_pulses*num_selected, 1);
for i=1:num_pulses
    cl_seq((i-1)*num_selected + 1:(i*num_selected), 1) = [modsig(i); single_pulse(2:end)]; 
end
md_seq = [mode; single_pulse(2:end)];
pw_seq = [phase_width; single_pulse(2:end)];
pg_seq = [phase_gap; single_pulse(2:end)];
per_seq = repmat(period, num_selected, 1);

ch_mag.channels   = repmat(ch_seq, num_pulses, 1);
ch_mag.magnitudes = cl_seq;
ch_mag.modes        = repmat(md_seq, num_pulses, 1);
ch_mag.phase_widths = repmat(pw_seq, num_pulses, 1);
ch_mag.phase_gaps   = repmat(pg_seq, num_pulses, 1);
ch_mag.periods      = repmat(per_seq, num_pulses, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
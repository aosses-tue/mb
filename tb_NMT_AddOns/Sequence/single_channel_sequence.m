function seq = single_channel_sequence( el, stim_rate, fm, range, params)
% function seq = single_channel_sequence( el, stim_rate, fm, range, params)
%
% el            - electrode number
% stim_rate     - stimulation rate [pps]
% range         - [min max] electrode to be displayed
% params        - structure with parameters for new sequence (optional)
%
% Comments by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    return;
end

if nargin < 5
    params = [];
end

params  = Ensure_field(params,'duration'   , 1); % in seconds
params  = Ensure_field(params,'phase_width',25); % in microseconds
params  = Ensure_field(params,'phase_gap'  , 8); % in microseconds 
params  = Ensure_field(params,'constant_current_level',0.8); % in percentage of DR
params  = Ensure_field(params,'mode', 103); % 103 = MP1+2
params  = Ensure_field(params,'NumMaxima',8); 

dur     = params.duration; % 1 second
dr      = params.constant_current_level;
mode    = params.mode;
pw      = params.phase_width;
pg      = params.phase_gap;
nb      = params.NumMaxima;

cm = Gen_single_channel_ch_mag_am_sequence(dur, el, dr, stim_rate, mode, pw, pg, nb, fm);
p = [];
p.analysis_rate = stim_rate;
p.channel_stim_rate = stim_rate;
p.num_selected = nb;
p = Channel_mapping_proc(p);
seq = Channel_mapping_proc(p, cm);
Plot_sequence(seq, '', range);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
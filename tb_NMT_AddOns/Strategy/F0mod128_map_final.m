function p = F0mod128_map_final(p)

% ACE_map: Calculate ACE map parameters.
% To perform processing, use:
%   q = Process(p, audio)
%
% p_out = ACE_map(p_in)
%
% p_in:  A struct containing the clinical parameters.
%          Any fields omitted will be set to default values.
% p_out: A struct containing the clinical and derived parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Authors: Matthias Milczynski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 0
%%
    p = [];
%%
end;
%%
p = Ensure_field(p, 'map_name', 'F0mod128');
p = Ensure_field(p, 'channel_stim_rate', 1800);
p = Ensure_field(p, 'analysis_rate', 1800);
p = Ensure_field(p, 'num_selected', 8);
p = Ensure_field(p, 'comfort_levels', ones(22, 1)*255);
p = Ensure_field(p, 'weight_strat', 'PS');
p = Ensure_field(p, 'block_length', 128);
p = Ensure_field(p, 'FAT', 'ACE');
p = Ensure_field(p, 'penv_filter_ord', 32); 
p = Ensure_field(p, 'penv_filter_cof', 60); % cutoff frequency in Hz
p = Ensure_field(p, 'sr', 16e3);
p = Ensure_field(p, 'gains_voiced_dB', 0);
p = Ensure_field(p, 'gains_unvoiced_dB', 0);
p = Ensure_field(p, 'magnitude_scaling', 1);
p = Ensure_rate_params(p);
p = Append_front_end_processes(p);

%%
p = Append_process(p, 'Autoc_opt2_Alt_proc');
%%
p = Append_process(p, 'FFT_filterbank_cell_in_proc');
%%
p = Append_process(p, 'Power_sum_envelope_cell_in_proc');
%%
p = Append_process(p, 'LPF_alt_order_proc');
%%
p = Append_process(p, 'AM_alt_order_proc');
%%
p = Append_process(p, 'Gain_voiced_alt_order_proc');
%%
p = Append_process(p, 'Gain_unvoiced_alt_order_proc');
%%
p = Append_process(p, 'Reject_smallest_proc');
p = Append_process(p, 'LGF_proc');

if (p.channel_stim_rate ~= p.analysis_rate)
	p = Append_process(p, 'Resample_FTM_proc');
end;

p = Append_process(p, 'Collate_into_sequence_proc');
p = Append_process(p, 'Scale_magnitude_proc');
p = Append_process(p, 'Channel_mapping_proc');
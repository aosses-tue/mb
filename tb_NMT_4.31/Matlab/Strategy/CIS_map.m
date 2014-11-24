function p = CIS_map(p)

% CIS_map: Calculate CIS map parameters.
% To perform processing, use:
%   q = Process(p, audio)
%
% p_out = CIS_map(p_in)
%
% p_in:  A struct containing the clinical parameters.
%          Any fields omitted will be set to default values.
% p_out: A struct containing the clinical and derived parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
	p = [];
end

p = Ensure_field(p, 'map_name', 'CIS');
p = Ensure_field(p, 'num_bands', 12);
p.num_selected = p.num_bands;

p = Ensure_rate_params(p);
p = Append_front_end_processes(p);
p = Append_process(p, 'FFT_filterbank_proc');
p = Append_process(p, 'Power_sum_envelope_proc');
p = Append_process(p, 'Gain_proc');
p = Append_process(p, 'LGF_proc');

if (p.channel_stim_rate ~= p.analysis_rate)
	p = Append_process(p, 'Resample_FTM_proc');
end

p = Append_process(p, 'Collate_into_sequence_proc');
p = Append_process(p, 'Channel_mapping_proc');

function p = ACE_map(p)
% function p = ACE_map(p)
%
% Calculate ACE map parameters.
% To perform processing, use:
%   q = Process(p, audio)
%
% p_out = ACE_map(p_in)
%
% p_in:  A struct containing the clinical parameters.
%          Any fields omitted will be set to default values.
% p_out: A struct containing the clinical and derived parameters.
%
% Edited by Alejandro Osses, ExpORL (added extra default parameters for
% compatibility with F0mod NMT Add-ons
%
% Default parameters:           value
%   -----------------------------------------------------------------------
% 	map_name                    'ACE'
%   scale_LGF_input (ExpORL)    0       if 1, an additional scaling is applied to the Input to the LGF, in the same way it is done on the xPC system
%
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

if isfield(p,'processes')
    disp([mfilename '.m processes will be deleted and reassigned'])
    disp('press any button to continue')
    pause()
    p = rmfield(p,'processes');
end

p = Ensure_field(p, 'map_name', 'ACE');
p = Ensure_field(p, 'scale_LGF_input',0);

p = Ensure_rate_params(p);
p = Append_front_end_processes(p);
p = Append_process(p, 'FFT_filterbank_proc');
p = Append_process(p, 'Power_sum_envelope_proc');
p = Append_process(p, 'Gain_proc');
p = Append_process(p, 'Reject_smallest_proc');
p = Append_process(p, 'LGF_proc');

if (p.channel_stim_rate ~= p.analysis_rate)
	p = Append_process(p, 'Resample_FTM_proc');
end

p = Append_process(p, 'Collate_into_sequence_proc');
p = Append_process(p, 'Channel_mapping_proc');

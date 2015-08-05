function p = ACE_map_new(p)
% function p = ACE_map_new(p)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%     $Archive: /Nucleus_MATLAB_Toolbox/Public/Strategy/ACE_map.m $
%    $Revision: 11 $
%        $Date: 17/01/05 2:49p $
%      Authors: Brett Swanson, Matthias Milczynski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 0
    p = [];
end;

p = Ensure_field(p, 'map_name'          , 'ACE');
p = Ensure_field(p, 'scale_LGF_input'   , 0); % Added by Alejandro
p = Ensure_field(p, 'audio_sample_rate' , 16000);
p = Ensure_field(p, 'channel_stim_rate' , 1800);
p = Ensure_field(p, 'analysis_rate'     , 1800);
p = Ensure_field(p, 'num_selected'      , 8);
p = Ensure_field(p, 'comfort_levels'    , ones(22, 1)*255);
p = Ensure_field(p, 'magnitude_scaling' , 1);

p = Ensure_field(p, 'source'            , 'Freedom'); % CP810 is not fully validated

p = Ensure_rate_params(p);
if p.processPreFiltered == 0
    p = Append_process(p,'Wav_proc');
    if p.normalize
        p = Append_process(p, 'Normalize_wav_proc');
    end
    if p.micfront
       p = Append_process(p,'FreedomMicResp_proc'); 
    end
else
     p = Append_process(p, 'Load_Prefiltered_proc');
end


switch p.filterbank
    case 'Freedom'
        p = Append_process(p, 'FFT_filterbank_proc');
        p = Append_process(p, 'Power_sum_envelope_proc');
    case '3G'
        p.crossover_freqs = ESPrit_crossover_frequencies(7, p.audio_sample_rate/2);
        p = Append_process(p, 'FFT_filterbank_proc');
        p = Append_process(p, 'Vector_sum_proc');
        p = Append_process(p,'Abs_proc');
    otherwise
        error('p.FILTERBANK can only be Freedom or 3G');
end

p = Append_process(p, 'Gain_proc');
p = Append_process(p, 'Reject_smallest_proc');
p = Append_process(p, 'LGF_proc');

if (p.channel_stim_rate ~= p.analysis_rate)
	p = Append_process(p, 'Resample_FTM_proc');
end;

p = Append_process(p, 'Collate_into_sequence_proc');
if ~p.QIC
    p = Append_process(p, 'Channel_mapping_proc');
end
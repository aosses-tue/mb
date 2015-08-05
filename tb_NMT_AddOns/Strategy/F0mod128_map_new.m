% function p = F0mod128_LPF_zero_ph_map(p)
%
% Latest F0mod128 version with zero-phase LPF  
%
% INPUTS
%   p   : map structure (can be also empty i.e. [])
%
% OUTPUTS
%   p   : Modified/Created map
% 
%      Author: Matthias Milczynski 2008

function p = F0mod128_map_new(p)

if nargin == 0
    p = [];
end;

%%
p = Ensure_field(p, 'rms_gain_dB', 0);
p = Ensure_field(p, 'map_name', 'F0m');
p = Ensure_field(p, 'channel_stim_rate', 1800);
p = Ensure_field(p, 'analysis_rate', 1800);
p = Ensure_field(p, 'num_selected', 8);
p = Ensure_field(p, 'comfort_levels', ones(22, 1)*255);
p = Ensure_field(p, 'weight_strat', 'PS');
p = Ensure_field(p, 'block_length', 128);
p = Ensure_field(p, 'FAT', 'ACE');
p = Ensure_field(p, 'penv_filter_ord', 8); 
p = Ensure_field(p, 'penv_filter_cof', 60); % cutoff frequency in Hz
p = Ensure_field(p, 'sr', 16e3);
p = Ensure_field(p, 'gains_voiced_dB', 0);
p = Ensure_field(p, 'gains_unvoiced_dB', 0);
p = Ensure_field(p, 'magnitude_scaling', 1);
p = Ensure_field(p, 'plot_modulator', 0);
% Can be sawtooth as well (according to Green et al. 2004, 2005)
p = Ensure_field(p, 'modshape', 'sine'); 
p = Ensure_field(p, 'adjust_absolute_level', 0);
p = Ensure_field(p, 'f0File', 0);
p = Ensure_field(p, 'processPreFiltered',0);
p = Ensure_field(p, 'normalize',0);
p = Ensure_field(p, 'micfront',0);
p = Ensure_field(p, 'filterbank','freedom');

p = LGF_proc(p);

if p.processPreFiltered == 0
    p = Append_process(p,'Wav_proc');
    p = Append_process(p, 'Adjust_speech_rms_proc');
    if p.normalize
        p = Append_process(p, 'Normalize_wav_proc');
    end
    if p.micfront
       p = Append_process(p,'FreedomMicResp_proc'); 
    end
else
     p = Append_process(p, 'Load_Prefiltered_proc');
end

if (p.f0File == 0)
    p = Append_process(p, 'Autoc_opt2_Alt_proc');
else
    p = Append_process(p, 'Read_F0_from_file_proc');
end

switch p.filterbank
    case 'freedom'
        p = Append_process(p, 'FFT_filterbank_cell_in_proc');
        p = Append_process(p, 'Power_sum_envelope_cell_in_proc');
    case '3G'
        p.audio_sample_rate = 16000;
        p.crossover_freqs = ESPrit_crossover_frequencies(7, p.audio_sample_rate/2);
        p = Append_process(p, 'FFT_filterbank_cell_in_proc');
        p = Append_process(p, 'Vector_sum_cell_in_proc');
        p = Append_process(p,'Abs_cell_in_proc');
    otherwise
        error('p.filterbank can only be freedom or 3G');
end
p = Append_process(p, 'LPF_zero_ph_proc');
p = Append_process(p, 'Gain_voiced_alt_order_proc');
p = Append_process(p, 'Gain_unvoiced_alt_order_proc');
p = Append_process(p, 'AM_BS_corr_proc');
p = Append_process(p, 'Reject_smallest_proc');
p = Append_process(p, 'LGF_proc');
if (p.channel_stim_rate ~= p.analysis_rate)
	p = Append_process(p, 'Resample_FTM_proc');
end
p = Append_process(p, 'Collate_into_sequence_proc');
p = Append_process(p, 'Scale_magnitude_proc');
if ~p.qic
    p = Append_process(p, 'Channel_mapping_proc');
end
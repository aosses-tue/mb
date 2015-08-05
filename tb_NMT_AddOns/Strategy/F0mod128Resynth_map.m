function p = F0mod128Resynth_map(p)
% function p = F0mod128Resynth_map(p)
%
% Latest F0mod128 version with zero-phase LPF  
%
% INPUTS
%   p   : map structure (can be also empty i.e. [])
%
% OUTPUTS
%   p   : Modified/Created map
%
% Dependencies:
%   setupStandardParamsF0Resynth
%
% Author: Matthias Milczynski 2009, edited by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    p = [];
end;

try
    p = setupStandardParamsF0Resynth(p);
catch
    addpath('/home/alejandro/Documenten/MATLAB/MATLAB_svn/NMTAddOns/Utility/')
    p = setupStandardParamsF0Resynth(p);
end
%%
p = Ensure_field(p, 'map_name', 'F0mod128');
p = Ensure_field(p, 'block_length', 128);
p = Ensure_field(p, 'penv_filter_ord', 8); 
p = Ensure_field(p, 'penv_filter_cof', 60); % cutoff frequency in Hz
p = Ensure_field(p, 'sr', 16e3);
p = Ensure_field(p, 'gains_voiced_dB', 0);
p = Ensure_field(p, 'gains_unvoiced_dB', 0);
p = Ensure_field(p, 'magnitude_scaling', 1);
p = Ensure_field(p, 'plot_modulator', 0);
% Can be sawtooth as well (according to Green et al. 2004, 2005)
p = Ensure_field(p, 'modshape', 'sine'); 
p = LGF_proc(p);

if p.processPreFiltered == 0
    %p = Append_process(p,'Wav_proc');
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
    case 'Freedom'
        p = Append_process(p, 'FFT_filterbank_cell_in_proc');
        p = Append_process(p, 'Power_sum_envelope_cell_in_proc');
    case '3G'
        p.audio_sample_rate = 16000;
        p.crossover_freqs = ESPrit_crossover_frequencies(7, p.audio_sample_rate/2);
        p = Append_process(p, 'FFT_filterbank_cell_in_proc');
        p = Append_process(p, 'Vector_sum_cell_in_proc');
        p = Append_process(p,'Abs_cell_in_proc');
    otherwise
        error('p.FILTERBANK can only be Freedom or 3G');
end
p = Append_process(p, 'LPF_zero_ph_proc');
p = Append_process(p, 'Gain_voiced_alt_order_proc');
p = Append_process(p, 'Gain_unvoiced_alt_order_proc');
p = Append_process(p, 'AM_BS_corr_proc');
p = Append_process(p, 'Reject_smallest_proc');
p = Append_process(p, 'ExpORLResynth_proc');
p = Append_process(p, 'Cos_ramp_proc');
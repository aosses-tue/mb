function p = F0mod128_map_new(p)

% function p = F0mod128_map_new(p)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:
%       Wav_proc (NMT), Read wav file and resample if necessary
%       Reject_smallest_proc (NMT)
%       LGF_proc (NMT), Loudness Growth.
%       Adjust_speech_rms_proc.m
%       Normalize_wav_proc.m
%       FreedomMicResp_proc.m
%       Load_Prefiltered_proc.m
%       Read_F0_from_file_proc.m
%       Autoc_opt2_Alt_proc
%       ESPrit_crossover_frequencies (NMT).
%       Vector_sum_cell_in_proc
%       Abs_cell_in_proc
%       LPF_zero_ph_proc                  ../MM/src/Matlab/DSP/NMTAddOns/FTM/
%       Gain_voiced_alt_order_proc
%       Gain_unvoiced_alt_order_proc
%       AM_BS_corr_proc
%       Resample_FTM_proc (NMT)
%       Collate_into_sequence_proc (NMT)
%       Scale_magnitude_proc              ../MM/src/Matlab/DSP/NMTAddOns/Sequence/
%       Channel_mapping_proc (NMT)
%       Gain_cell_in_proc.m               ../MM/src/Matlab/DSP/NMTAddOns/Processing/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % Example, loads default p-structure for F0mod:
%   addpath('/home/alejandro/Documenten/MM/src/Matlab/DSP/NMTAddOns/Processing/')  
%   addpath('/home/alejandro/Documenten/MM/src/Matlab/DSP/NMTAddOns/FTM/')
%   addpath('/home/alejandro/Documenten/MM/src/Matlab/DSP/NMTAddOns/Sequence/')
%   p = F0mod128_map_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    p = [];
end;

if isfield(p,'test')
    return;
end

%%
p = Ensure_field(p, 'input_scaling'     ,1);%
p = Ensure_field(p, 'processPreFiltered', 0);%
p = Ensure_field(p, 'gains_dB'          , 0);%
p = Ensure_field(p, 'filterbank'        ,'freedom');%
p = Ensure_field(p, 'micfront'          ,1);% % by default we will add SP12 Mic response

p = Ensure_field(p, 'rms_gain_dB'       , 0);
p = Ensure_field(p, 'map_name'          , 'F0m');
p = Ensure_field(p, 'channel_stim_rate' , 1800);
p = Ensure_field(p, 'analysis_rate'     , 1800); 
p = Ensure_field(p, 'num_selected'      , 8);
p = Ensure_field(p, 'comfort_levels'    , ones(22, 1)*255);
p = Ensure_field(p, 'weight_strat'      , 'PS');
p = Ensure_field(p, 'block_length'      , 128);
p = Ensure_field(p, 'FAT'               , 'ACE');
p = Ensure_field(p, 'penv_filter_ord'   , 8); 
p = Ensure_field(p, 'penv_filter_cof'   , 60); % cutoff frequency in Hz
p = Ensure_field(p, 'sr'                , 16e3);
p = Ensure_field(p, 'gains_voiced_dB'   , 0);
p = Ensure_field(p, 'gains_unvoiced_dB' , 0);
p = Ensure_field(p, 'magnitude_scaling' , 1);
p = Ensure_field(p, 'plot_modulator'    , 0);
% Can be sawtooth as well (according to Green et al. 2004, 2005)
p = Ensure_field(p, 'modshape'          , 'sine'); 
p = Ensure_field(p, 'adjust_absolute_level', 0);
p = Ensure_field(p, 'f0File'            , 0);

p = LGF_proc(p);

if p.processPreFiltered == 0 % Equivalent to 'Append_front_end_processes'
    
    p = Ensure_field(p,'micfront' ,1); % by default we will add SP12 Mic response
    if isfield(p,'processes')
        p = rmfield(p,'processes');
    end
    if p.micfront
        p = Append_front_end_processes(p);    
    end
    % p = Append_process(p, 'Adjust_speech_rms_proc'); % Calibration to p.rms_gain_dB (default = 0 dB)
    
else
     p = Append_process(p, 'Load_Prefiltered_proc'); % AOV
end
    
if (p.f0File == 0)
    p = Append_process(p, 'Autoc_opt2_Alt_proc');
else
    p = Append_process(p, 'Read_F0_from_file_proc'); % If preprocessing
end

p = Append_process(p, 'FFT_filterbank_cell_in_proc');

switch p.filterbank
    case 'freedom'
        p = Append_process(p, 'Power_sum_envelope_cell_in_proc');
    case '3G'
        p.audio_sample_rate = 16000;
        p.crossover_freqs = ESPrit_crossover_frequencies(7, p.audio_sample_rate/2);
        p = Append_process(p, 'Vector_sum_cell_in_proc');
        p = Append_process(p,'Abs_cell_in_proc');
    otherwise
        error('p.filterbank can only be freedom (includes CP810) or 3G');
end

if p.gains_dB == 1
    p = Append_process(p, 'Gain_cell_in_proc'); % Applying 0 dB, but just to have the same procedures than the default ACE
end
   
p = Append_process(p, 'LPF_zero_ph_proc');
p = Append_process(p, 'AM_BS_corr_proc');
p = Append_process(p, 'Reject_smallest_proc');
p = Append_process(p, 'LGF_proc');

if (p.channel_stim_rate ~= p.analysis_rate)
	p = Append_process(p, 'Resample_FTM_proc');
end

p = Append_process(p, 'Collate_into_sequence_proc');
p = Append_process(p, 'Scale_magnitude_proc');

if ~isfield(p,'qic')
    p = Append_process(p, 'Channel_mapping_proc');
end   

function p = F0mod_map_201401(p)
% function p = F0mod128_map_201401(p)
%
% Latest F0mod128 version with zero-phase LPF  
%
% INPUTS
%   p   : map structure (can be also empty i.e. [])
%
% OUTPUTS
%   p   : Modified/Created map
% 
% % Example, loads default p-structure for F0mod:
%   addpath('/home/alejandro/Documenten/MM/src/Matlab/DSP/NMTAddOns/Processing/')  
%   addpath('/home/alejandro/Documenten/MM/src/Matlab/DSP/NMTAddOns/FTM/')
%   addpath('/home/alejandro/Documenten/MM/src/Matlab/DSP/NMTAddOns/Sequence/')
%   p = F0mod_map_201401;
%
% Author: Matthias Milczynski 2008
% Edited by Alejandro Osses, ExpORL, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    p = [];
end;

if isfield(p,'test')
    return;
end

%%
p = Ensure_field(p, 'map_name'          , 'F0m');
p = Ensure_field(p, 'input_scaling'     ,1);%
p = Ensure_field(p, 'processPreFiltered', 0);%
p = Ensure_field(p, 'gains_dB'          , 0);%
p = Ensure_field(p, 'filterbank'        ,'Freedom');
p = Ensure_field(p, 'micfront'          ,1);% % by default we will add SP12 Mic response

p = Ensure_field(p, 'rms_gain_dB'       , 0);
p = Ensure_field(p, 'channel_stim_rate' , 900);
p = Ensure_field(p, 'analysis_rate'     , 900); 
p = Ensure_field(p, 'num_selected'      , 8);
p = Ensure_field(p, 'comfort_levels'    , ones(22, 1)*255);
p = Ensure_field(p, 'weight_strat'      , 'PS');
p = Ensure_field(p, 'block_length'      , 128);

p = Ensure_field(p, 'autoc_block_size'  , 352);
p = Ensure_field(p, 'autoc_block_shift' , 144);

p = Ensure_field(p, 'FAT'               , 'ACE');
p = Ensure_field(p, 'penv_filter_ord'   , 8); 
p = Ensure_field(p, 'penv_filter_cof'   , 60); % cutoff frequency in Hz
% p = Ensure_field(p, 'sr'                , 15659);
p = Ensure_field(p, 'audio_sample_rate' , 15659);
p = Ensure_field(p, 'gains_voiced_dB'   , 0);
p = Ensure_field(p, 'gains_unvoiced_dB' , 0);
p = Ensure_field(p, 'magnitude_scaling' , 1);
p = Ensure_field(p, 'plot_modulator'    , 0);
% Can be sawtooth as well (according to Green et al. 2004, 2005)
p = Ensure_field(p, 'modshape'          , 'sine'); 
p = Ensure_field(p, 'adjust_absolute_level', 0);
p = Ensure_field(p, 'f0File'            , 0);
p = Ensure_field(p, 'silence_thr'       , 0.01);

p = LGF_proc(p);

if p.processPreFiltered == 0 % Equivalent to 'Append_front_end_processes'
    
    p = Ensure_field(p,'micfront' ,1); % by default we will add SP12 Mic response
    if isfield(p,'processes')
        p = rmfield(p,'processes');
    end
        
else
     p = Append_process(p, 'Load_Prefiltered_proc'); % AOV
end
    
if (p.f0File == 0)
    p = Append_process(p, 'Autoc_opt2_Alt_proc');
else
    p = Append_process(p, 'Read_F0_from_file_proc'); % If preprocessing
end

p = Append_process(p, 'FFT_filterbank_cell_in_proc');
p = Append_process(p, 'Power_sum_envelope_cell_in_proc');


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

p = Append_process(p, 'Channel_mapping_proc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

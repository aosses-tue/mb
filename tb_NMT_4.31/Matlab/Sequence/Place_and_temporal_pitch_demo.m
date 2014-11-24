% Place_and_temporal_pitch_demo: Plot sequences that vary in place or temporal pitch.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Cochlear Ltd
%   $Header: //dsp/Nucleus/_Releases/NMT_4.30/Matlab/Sequence/Place_and_temporal_pitch_demo.m#1 $
% $DateTime: 2008/03/04 14:27:13 $
%   $Change: 86418 $
%   Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [];
x.stimulus_duration	= 0.1;

p = [];
p.num_bands         = 22;
p.num_selected		= 6;
p.channel_stim_rate = 2400;
p.implant_stim_rate = 2400;
p = Append_process(p, 'Channel_mapping_proc');

x_e = Gen_channel_rank_sequences(p, x);
x_e = Process_stimuli(p, x_e);
Plot_sequence(x_e.stimuli, 'Place pitch');

x_spp = x;
x_spp = Gen_SPP_rank_sequences(p, x_spp);
x_spp = Process_stimuli(p, x_spp);
Plot_sequence(x_spp.stimuli, 'Rate pitch', 1:22);

x_am = x;
x_am = Gen_AM_rank_sequences(p, x_am);
x_am = Process_stimuli(p, x_am);
Plot_sequence(x_am.stimuli, 'Modulation pitch', 1:22);

x_mpp = x;
x_mpp = Gen_HWG_sync_rank_sequences(p, x_mpp);
x_mpp = Process_stimuli(p, x_mpp);
Plot_sequence(x_mpp.stimuli, 'Modulation (MPP) pitch', 1:22);


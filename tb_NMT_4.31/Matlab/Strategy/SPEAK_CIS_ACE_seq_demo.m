% SPEAK_CIS_ACE_seq_demo: Plot sequences of "choice" for each strategy.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Cochlear Ltd
%   $Header: //dsp/Nucleus/_Releases/NMT_4.30/Matlab/Strategy/SPEAK_CIS_ACE_seq_demo.m#1 $
% $DateTime: 2008/03/04 14:27:13 $
%   $Change: 86418 $
%   Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

audio = 'choice';

p = [];
p.source = 'HS8';
p.num_bands = 22;
p.audio_sample_rate = 16000;

p_spk = p;
p_spk.num_selected = 6;
p_spk.channel_stim_rate = 250;
p_spk = ACE_map(p_spk);

p_ace = p;
p_ace.num_selected = 12;
p_ace.channel_stim_rate = 1200;
p_ace = ACE_map(p_ace);

p_cis = p;
p_cis.num_bands = 6;
p_cis.electrodes = (22:-4:1)';
p_cis.channel_stim_rate = 2400;
p_cis = CIS_map(p_cis);

q_spk = Process(p_spk, audio);
q_cis = Process(p_cis, audio);
q_ace = Process(p_ace, audio);

Plot_sequence({q_cis, q_spk, q_ace}, {'CIS', 'SPEAK', 'ACE'}, 1:22);

% For pasting figures into Word doc:
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [1, 1, 16, 10]);

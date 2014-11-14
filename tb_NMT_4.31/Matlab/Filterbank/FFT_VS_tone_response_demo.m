% FFT_VS_tone_response_demo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ps.crossover_freqs = [
 142
 297
 480
 660
 840
1018
1196
1375
1571
1805
2073
2382
2736
3175
3720
4358
5101
5982
7008
];

ps.buffer_opt = 'nodelay';	% Skip initial transient
ps = Append_process(ps, 'FFT_filterbank_proc'); 
ps = Append_process(ps, 'Vector_sum_proc'); 
ps = Append_process(ps, 'Abs_proc'); 

GUI_FTM(ps, ps.freq_response);
Plot_freq_response(ps.response_freqs, ps.freq_response, 'log');

Fs = ps.bin_freq/8;
freq_vec = 0:Fs:7500;
[output, response, max_gain] = Tone_response(ps, freq_vec, 2*ps.block_length);

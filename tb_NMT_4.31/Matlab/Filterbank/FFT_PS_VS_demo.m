% FFT_PS_VS_demo: Comparison of power-sum and vector-sum envelopes for FFT filterbank.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filterbank parameters (defaults):
p = [];
p = Append_front_end_processes(p);
p = Append_process(p, 'FFT_filterbank_proc');

p_ps = p;
p_ps = Append_process(p_ps, 'Power_sum_envelope_proc'); 

p_vs = p;
p_vs = Append_process(p_vs, 'Vector_sum_proc'); 
p_vs = Append_process(p_vs, 'Abs_proc');

% Process and display:
if ~exist('audio', 'var')
	audio = 'asa';
end
u_ps = Process(p_ps, audio);
u_vs = Process(p_vs, audio);

GUI_FTM(p_ps, {u_ps, u_vs}, {'Env power sum', 'Env vector sum'});

p_tail = p;
p_tail.processes = {};
p_tail = Append_process(p_tail, 'Reject_smallest_proc');
p_tail = Append_process(p_tail, 'LGF_proc');
p_tail = Append_process(p_tail, 'Collate_into_sequence_proc');

q_ps = Process(p_tail, u_ps);
q_vs = Process(p_tail, u_vs);
Plot_sequence({q_ps, q_vs}, {'Seq power sum', 'Seq vector sum'});

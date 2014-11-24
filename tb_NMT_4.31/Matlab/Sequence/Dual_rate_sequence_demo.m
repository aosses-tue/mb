% Dual_rate_sequence_demo: Demonstrates a sequence with different rates on two channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_pulses		= 1000;
carrier_rate	= 10000;

rq.channels		= repmat(0,  num_pulses, 1);
rq.magnitudes	= repmat(-1, num_pulses, 1);
rq.periods		= 1e6/carrier_rate; % microseconds

rate1 = 500;
N1 = round(carrier_rate/rate1);
k1 = 1:N1:num_pulses;
rq.channels  (k1) = 21;
rq.magnitudes(k1) = 1;

rate2 = 1000;
N2 = round(carrier_rate/rate2);
k2 = 2:N2:num_pulses;
rq.channels  (k2) = 22;
rq.magnitudes(k2) = 1;

p  = Channel_mapping_proc;
pq = Channel_mapping_proc(p, rq);

Plot_sequence(pq);

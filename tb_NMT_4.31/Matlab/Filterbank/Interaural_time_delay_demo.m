% Interaural_time_delay_demo: Shows a typical delay between left & right ears on a speech sample.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = wavread('asa');
delay = 0.3; % milliseconds
region = [138, 168];

n = 1:length(a);
t = (n-1)/16;	% milliseconds
figure;
hold on;
plot(t, a);
td = t - delay;
plot(td, a, 'r')
xlabel('Time (milliseconds)')
set(gca, 'XLim', region)
Window_title('Waveforms')
zoom xon

p = [];
p.audio_sample_rate		= 16000;
p.analysis_rate			= p.audio_sample_rate;
p.num_bands				=    22;

p = Append_process(p, 'FFT_filterbank_proc');
p = Append_process(p, 'Vector_sum_proc');

fb_out = Process(p, a);
%GUI_FTM(p, fb_out);

env = abs(fb_out);
channel = 4;
v = env(channel, :);
n = 1:length(v);
t = (n-1)/16;	% milliseconds
td = t - delay;
figure;
hold on;
plot(t, v);
plot(td, v, 'r');
xlabel('Time (milliseconds)')
set(gca, 'XLim', region)
Window_title('Filter envelopes')
zoom xon;


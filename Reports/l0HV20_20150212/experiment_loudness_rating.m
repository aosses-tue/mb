clear all
close all
% seed = 10;
% rng(seed); % run random number generator with seed.
sf = 44100;
dt = 1/sf;
f0 = 1000;
burst_duration = 0.25;

participant = 'Dik';
% participant = input('Give participant''s name please: ', 's');

fade_duration = 0.01;
inter_burst_interval = 0.05;
inter_pair_interval = 1.00;

dB_diff = -20:5:20; % array of differences in dB SPL
                    % between bursts of pair
dB_order = dB_diff(randperm(length(dB_diff)));

t_burst = 0:dt:burst_duration;
t_itp = 0:dt:inter_burst_interval;
t_ipp = 0:dt:inter_pair_interval;
inter_burst_pause = zeros(size(t_itp));
inter_pair_pause = zeros(size(t_ipp));

burst = 0.2*randn(1, length(t_burst));
burst = Fade(burst, sf, fade_duration);

tim = sprintf('%s: %4d %2d %2d; %2dh %2dm %2.0fs', participant, clock);
output_filename = ...
    sprintf('loudness_rating_%s %s.txt', participant, tim);
for j = 1:length(output_filename)
    if output_filename(j) == ' '
        output_filename(j) = '_';
    end
end
fp = fopen(output_filename, 'wt');
fprintf(fp, '%s\n', tim);
fprintf('%s\n', tim);
stim = [];
stim_text = [];
for k = 1:length(dB_order)
    inc = 10^(dB_order(k)/20);
    txt = sprintf('%7.0f', dB_order(k));
    stim_text = [stim_text txt]; %#ok<*AGROW>
    second_burst = inc*burst;
    stim = [stim burst inter_burst_pause second_burst inter_pair_pause]; 
    fprintf(fp, '%4.0f', dB_order(k));
    fprintf('%4.0f', dB_order(k));
end
fprintf(fp, '  dB\n');
fprintf('  dB\n');
sound(stim, sf)
fclose(fp);
disp('Thank you for your attention')

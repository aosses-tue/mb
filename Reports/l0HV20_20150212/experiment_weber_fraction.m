clear all
close all
% seed = 10;
% rng(seed); % run random number generator with seed.
sf = 44100;
dt = 1/sf;
f0 = 1000;
tone_duration = 0.25;

participant = input('Give your name please: ', 's');

fade_duration = 0.01; %sec
inter_tone_interval = 0.25; %sec
inter_pair_interval = 2.00; %sec
inter_set_interval = 3; %sec

nr_of_pairs = 10; % nr of pairs presented in one run
dB_range = -20:10:20; % range of intensities
dB_diff = [5:-1:0]; % decreasing array of differences in dB SPL
                    % between tones of pair
dB_order = dB_range(randperm(length(dB_range)));
a_max = 10^(-max(dB_diff)/20);
a_range = a_max * 10.^((dB_order-max(dB_order))/20);

t_tone = 0:dt:tone_duration;
t_itp = 0:dt:inter_tone_interval;
t_ipp = 0:dt:inter_pair_interval;
t_irp = 0:dt:inter_set_interval;
inter_tone_pause = zeros(size(t_itp));
inter_pair_pause = zeros(size(t_ipp));
inter_run_pause = zeros(size(t_irp));

tone = sin(2*pi*f0*t_tone);
tone = Fade(tone, sf, fade_duration);

tim = sprintf('%4d %2d %2d; %2dh %2dm %2.0fs', clock);
output_filename = ...
    sprintf('weber_fraction_%s %s.txt', participant, tim);
for j = 1:length(output_filename)
    if output_filename(j) == ' '
        output_filename(j) = '_';
    end
end
fp = fopen(output_filename, 'wt');

fprintf(fp, '%s', tim);
fprintf('%s', tim);
for k = 1:length(a_range)
    first_tone = a_range(k) * tone;
    stim = [];
    stim_text = [];
    for i = 1:length(dB_diff)
        inc = 10^(dB_diff(i)/20);
        up_down = round(rand(1, nr_of_pairs));
        up_down_text = sprintf('%4d', up_down);
        up_down_text = strrep(up_down_text, '1', 'd');
        up_down_text = strrep(up_down_text, '0', 'u');
        txt = sprintf('%7.2f %s\n', dB_diff(i), up_down_text);
        if i == 1
            stim_text = txt;
        else
            stim_text = [stim_text txt];
        end
        second_tone = inc*first_tone;
        run = [];
        for j = 1:length(up_down)
            if up_down(j) == 0
                tone_pair = ...
                    [first_tone inter_tone_pause second_tone];
            else
                tone_pair = ...
                    [second_tone inter_tone_pause first_tone];
            end
            run = [run tone_pair inter_pair_pause]; %#ok<*AGROW>
        end
        stim = [stim run inter_run_pause];
    end
    fprintf(fp, '\nrelative level: %7.3f dB\n', dB_order(k));
    fprintf(fp, '%s\n', stim_text);
    fprintf('\nrelative level: %7.3f dB\n', dB_order(k));
    fprintf('%s\n', stim_text);
    sound(stim, sf)
    pause
end
fclose(fp);
disp('Last run')
disp('Thank you for your attention')

clear all
close all
% seed = 10;
% rng(seed); % set random number generator with seed.
sf = 44100;
dt = 1/sf;
f0 = 1000;
tone_duration = 0.25;

fade_duration = 0.01;
inter_tone_interval = 0.05;
inter_pair_interval = 0.1;
inter_set_interval = 0.2;

nr_of_pairs = 1;
dB_range = -20:10:20;
dB_range = dB_range(randperm(length(dB_range)));
dB_diff = 10;
a_max = 10^(-max(dB_diff)/20);
a_range = a_max * 10.^((dB_range-max(dB_range))/20);

t_tone = 0:dt:tone_duration;
t_itp = 0:dt:inter_tone_interval;
t_ipp = 0:dt:inter_pair_interval;
t_isp = 0:dt:inter_set_interval;
inter_tone_pause = zeros(size(t_itp));
inter_pair_pause = zeros(size(t_ipp));
inter_set_pause = zeros(size(t_isp));

fprintf('\n----------------------------------------------- \n')
fprintf('ADJUSTMENT OF SOUND LEVEL\n')
fprintf('\n----------------------------------------------- \n')
fprintf('You will hear tone pairs varying over the intensity \n')
fprintf('range within which they will be presented during \n')
fprintf('the experiment. The difference in intensity between \n')
fprintf('the two toens of a tone pair is now maximum.\n\n')
fprintf('Adjust the sound level until the softest tone pair  \n')
fprintf('is still audible, but the loudest pairs are not \n')
fprintf('too loud. \n')
fprintf('To repeat the tone pairs, press ''y'' or <enter>.\n')
fprintf('When you are ready, press ''n'' \n')
fprintf('\n')
c = input('Now press <enter>!  ', 's');
if isempty(c)
    c = 'y';
end

tone = sin(2*pi*f0*t_tone);
tone = Fade(tone, sf, fade_duration);
while c == 'y' || c == 'Y'
    stim = [];
    for k = 1:length(a_range)
        first_tone = a_range(k) * tone;
        for i = 1:length(dB_diff)
            inc = 10^(dB_diff(i)/20);
            up_down = round(rand(1, nr_of_pairs));
            % up_down_text = sprintf('%4d', up_down);
            % txt = sprintf('%7.2f %s', dB_diff(i), up_down_text);
            % if i == 1
            %     stim_text = txt;
            % else
            %     stim_text = [stim_text; txt];
            % end
            second_tone = inc*first_tone;
            set = [];
            for j = 1:length(up_down)
                if up_down(j) == 0
                    tone_pair = ...
                        [first_tone inter_tone_pause second_tone];
                else
                    tone_pair = ...
                        [second_tone inter_tone_pause first_tone];
                end
                set = [set tone_pair inter_pair_pause]; %#ok<*AGROW>
            end
            stim = [stim set inter_set_pause];
        end
        % fprintf('\ndB range: %7.3f\n', dB_range(k))
        % disp(stim_text)
    end
    sound(stim, sf)
    disp(' ')
    c = input('Repeat (y/n)? ', 's');
    if isempty(c)
        c = 'y';
    end
end
fprintf('\n')
fprintf('Do not change the sound level \n')
fprintf('during the experiment!\n\n')
fprintf('Now do the experiment. \n')


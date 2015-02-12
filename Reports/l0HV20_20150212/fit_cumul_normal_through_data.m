clear all
close all

% Ask for the input level of this data set
level = input('The level difference from average is (-20:10:20): ');

set(gcf, 'Position', [240 100 800 500], ...
    'PaperPosition', [0.6 0.6 16 10], 'InvertHardCopy', 'off')
set(gca, 'Fontsize', 14, 'FontWeight', 'bold', ...
    'box', 'on', 'NextPlot', 'add')

% The differences in dB between tones of a pair
% used in the constant stimuli experiment
intensity_diffs = 0:2:10;

% The next line 18 simulates a possible outcome of your experiment.
% Replace the data with the actual outcome of your experiment.
% *************************************************************** %
% The array YDATA must contain the percentages correct you obtained.
ydata = [53.2   46.1   65.3   81.8   92.7  103.1];
% *************************************************************** %

% Plot the points of the measured psychometric function
plot(intensity_diffs, ydata, 'xk', 'MarkerSize', 12, 'LineWidth', 2)

% Estimate the parameters of the cumulative, normal fit to the data.
[estimates, model] = Fitcumulnormal(intensity_diffs, ydata);
mean = estimates(1);
stdev = estimates(2);

% Show the estimated mean and standard deviation of the
% distribution
disp(' ')
fprintf('Estimated mean:     %4.2f\n', mean);
fprintf('Estimated st. dev.: %4.2f\n\n', stdev);

% Range of x values for which the fitting curve will
% be calculated and the corresponding fit.
xrange = (intensity_diffs(1)-0.9):0.01:(intensity_diffs(end)+0.9);
fit = 50+50*normcdf(xrange, mean, stdev);

% Plot the estimated cumulative normal distribution through the
% obtained estimates of the psychometric function
plot(xrange, fit, 'r', 'LineWidth', 2)
axis([xrange(1)-1 xrange(end)+1 41 109])
xlabel('intensity diffence (dB)')
ylabel('percentage correct (%)')
text(-1, 100, sprintf('level diff. from average: %4.0f dB', level), ...
    'Fontsize', 14, 'FontWeight', 'bold')
text(-1, 90, sprintf('mean:      %4.2f dB', mean), ...
    'Fontsize', 14, 'FontWeight', 'bold')
text(-1, 85, sprintf('st. dev.:   %4.2f dB', stdev), ...
    'Fontsize', 14, 'FontWeight', 'bold')

tim = sprintf('%4d %2d %2d; %2dh %2dm %2.0fs', clock);
if level == 0
    filename = ...
        sprintf('%s %s %s.jpg', mfilename, '__0_', tim);
else
    filename = ...
        sprintf('%s %s %s.jpg', mfilename, sprintf('%3.0d', level), tim);
end
for j = 1:length(filename)
    if filename(j) == ' '
        filename(j) = '_';
    end
end
print(gcf, '-djpeg', sprintf('%s', filename))


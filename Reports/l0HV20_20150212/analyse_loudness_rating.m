clear all
close all

% Resultaten Loudness Scaling - Sound Design 2000
results = [200 100  50 100  50 200 150 150  25 100  75 175 200 100 125  15 100  25 125 150
    200  75  50 100  80 300 120 160  40  90  70 250 300  80 150  60  90  30 130 200
    200 100  50 100  50 300 150 200  50 100  50 200 300 100 200  50 100  25 100 300
    200  95  30 100  80 240 110 120  70 100  80 220 280  95 115  40  85  30 110 200
    250  95  60  95  50 250 150 175  40  95  70 225 250  95 125  40  90  50  75 225];

levels = [15  -5 -20   0 -10  20   5  10 -15   0 -10  15  20  -5  10 -15  -5 -20   5  15];
sl = repmat(levels, size(results, 1), 1);
yl = log10([10 25 50 100 250 500 1000]);

set(gcf, 'Position', [240 100 800 800], ...
    'PaperPosition', [0.6 0.6 15 10], 'InvertHardCopy', 'off')
set(gca, 'linewidth', 2, 'fontsize', 14, 'box', 'on', ...
    'YTick', yl, 'YTickLabel', 10.^yl, 'NextPlot', 'add')

logres = log10(results);
p = round(levels/5+5);
xlevels(p) = levels;
% calculate and plot quartiles
quartiles = zeros(9, 3);
for k = 1:9
    zz = logres(:, p == k);
    quartiles(k, :) = prctile(zz(:), [25 50 75]);
end
dx = 1.0;
for k = 1:9
    xk = 5*k-25;
    plot([xk xk xk], quartiles(k, :), 'LineWidth', 2)
    xp = [xk-dx xk+dx];
    plot(xp, [quartiles(k, 1) quartiles(k, 1)], 'LineWidth', 3)
    plot(xp, [quartiles(k, 2) quartiles(k, 2)], 'LineWidth', 3)
    plot(xp, [quartiles(k, 3) quartiles(k, 3)], 'LineWidth', 3)
end
plot(sl, logres, 'xk', 'MarkerSize', 9, 'LineWidth', 2);
axis([-25 25 1 3]);
xlabel('sound level (dB)', 'fontsize', 14);
ylabel('loudness rating ', 'fontsize', 14);
title('loudness rating experiment', 'fontsize', 14)
% calculate and plot regression line
[a0, a1] = LinRegress(sl, logres);
disp([a0 a1]);
expo = sprintf('exponent =%6.3f', 10*a1);
text(6, 1.25, expo, 'fontsize', 14)
text(-20, 2.7, sprintf('N = %2d', size(results, 1)), 'fontsize', 14)
x1 = -30;
x2 = 30;
y1 = a1*x1+a0;
y2 = a1*x2+a0;
plot([x1 x2], [y1 y2], '-r', 'linewidth', 2);
print(gcf, '-djpeg', sprintf('%s_1', mfilename))

figure
set(gcf, 'Position', [300 100 600 800], ...
    'PaperPosition', [0.6 0.6 9 15], 'InvertHardCopy', 'off')

% plot results on a log-log scale
subplot(3,1,3)
set(gca, 'linewidth', 2, 'fontsize', 10, 'fontweight', 'bold', ...
    'box', 'on', 'NextPlot', 'add', ...
    'YTick', yl, 'YTickLabel', 10.^yl, 'TickDir', 'Out')
plot(xlevels', quartiles(:, 2), '-x', 'linewidth', 2)
for j=1:length(xlevels),
    plot([xlevels(j) xlevels(j)], [quartiles(j, 1) quartiles(j, 3)])
end
plot([x1 x2], [y1 y2], '--r', 'linewidth', 1);
axis([-22 22 1 log10(330)])
xlabel('relative sound level (dB SPL)', ...
    'fontsize', 10, 'fontweight', 'bold')
ylabel('rating (log %)', 'fontsize', 10, 'fontweight', 'bold')

% plot results with a logarithmic ordinate
subplot(3,1,2)
yl = 0:100:300;
set(gca, 'linewidth', 2, 'fontsize', 10, 'fontweight', 'bold', ...
    'box', 'on', 'NextPlot', 'add', ...
    'YTick', yl, 'YTickLabel', yl, 'TickDir', 'Out')
plot(xlevels', 10.^quartiles(:, 2), '-x', 'linewidth', 2)
for j=1:length(xlevels),
    plot([xlevels(j) xlevels(j)], [10^quartiles(j, 1) 10^quartiles(j, 3)])
end
axis([-22 22 0 330])
xlabel('relative sound level (dB SPL)', ...
    'fontsize', 10, 'fontweight', 'bold')
ylabel('rating (%)', 'fontsize', 10, 'fontweight', 'bold')

subplot(3,1,1)
% plot results with a logarithmic abscissa
set(gca, 'linewidth', 2, 'fontsize', 10, 'fontweight', 'bold', ...
    'box', 'on', 'NextPlot', 'add', ...
    'YTick', yl, 'YTickLabel', yl, 'TickDir', 'Out')
ylevels = 20*exp(xlevels/10);
plot(ylevels, 10.^quartiles(:, 2), '-x', 'linewidth', 2)
for j=1:length(ylevels),
    plot([ylevels(j) ylevels(j)], [10^quartiles(j, 1) 10^quartiles(j, 3)])
end
axis([0 150 0 330])
xlabel('relative intensity (~W/m^2)', ...
    'fontsize', 10, 'fontweight', 'bold')
ylabel('rating (%)', 'fontsize', 10, 'fontweight', 'bold')
print(gcf, '-djpeg', sprintf('%s_2', mfilename))

disp('')

function r20160504_update_analysing_results
% function r20160504_update_analysing_results
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       See also: r20160429_experiments_WAE.m
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 04/05/2016
% Last update on: 04/05/2016 
% Last use on   : 04/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

dirmain   = ['D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20160503-Acoustics-TUe' delim];
dirthis   = ['Figures-MATLAB'      delim];
dirres    = [dirmain 'Experiments' delim];
diraudio  = [dirmain 'Sounds'      delim];

outputdir = [dirmain dirthis];
Mkdir(outputdir);
Mkdir(diraudio);

FontSize = 14;
bFig7    = 1;

h = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotOpts = [];
plotOpts.I_Width = 12;
plotOpts.FontSize = FontSize;

if bFig7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dir = 'D:\Databases\dir04-Psychoacoustics\WAE\saves-my-results\';
    % files = {'save-AO-silence-order.xml','save-AO-silence-1.xml'};
    dir = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\01-Results-summary\';
    
    files       = { 'save-LB-silence.xml', 'save-XK-silence.xml', 'save-AOLab-silence.xml'};
    files_noise = {'save-LB-in-noise.xml','save-XK-in-noise.xml','save-AOLab-in-noise.xml'};
    
    score_pool = [];
    score_pool_noise = [];
    
    for i = 1:length(files)
        [scoref_pair_tmp] = quick_Triads_WAE([dir files{i}]);
        score_pool = [score_pool; transpose(scoref_pair_tmp(:,1))];
    end
            
    label_piano={'0','2t1','2t2','4','6'};
    for i = 1:size(scoref_pair_tmp,1)
        p1 = scoref_pair_tmp(i,2);
        p2 = scoref_pair_tmp(i,3);
        label_pair{i} = sprintf('%s/%s',label_piano{p1},label_piano{p2});
    end
    
    score_avg = median(score_pool);
    scoreU    = prctile(score_pool,95);
    scoreL    = prctile(score_pool, 5);
        
    [label_pair idx] = sort(label_pair);
    score_avg        = score_avg(idx);
    scoreL = scoreL(idx);
    scoreU = scoreU(idx);
    
    [xx idx2]       = sort(score_avg,'descend'); % second ordering criterion
    score_avg       = score_avg(idx2);
    scoreL = scoreL(idx2);
    scoreU = scoreU(idx2);
    label_pair      = label_pair(idx2);
    
    errorU    = scoreU - score_avg;
    errorL    = score_avg - scoreL; 
    
    offsetx = 0.1;
    figure;
    xdata = 1:length(score_avg);
    errorbar(xdata-offsetx,score_avg,errorL,errorU,'ro','LineWidth',2); grid on, hold on
    Xlabel('Piano pair',FontSize)
    Ylabel('Most similar pair, Score[\%]',FontSize)
    
    xlim([min(xdata)-0.5 max(xdata)+0.5]);
    ylim([0 103])
    set(gca,'XTick',xdata);
    h(end+1) = gcf;
    h(end) = Figure2paperfigureT(h(end),1,1,plotOpts);
    set(gca,'XTickLabel',label_pair);
    set(gca,'FontSize',FontSize);
    hold on;
    legend('silence')
    
    % name = 'triadic-three-participants';
    % Saveas(h(end),sprintf('%s%s',outputdir,name),'epsc')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(files_noise)
        [scoref_pair_tmp] = quick_Triads_WAE([dir files_noise{i}]);
        score_pool_noise = [score_pool_noise; transpose(scoref_pair_tmp(:,1))];
    end
            
    label_piano={'0','2t1','2t2','4','6'};
    for i = 1:size(scoref_pair_tmp,1)
        p1 = scoref_pair_tmp(i,2);
        p2 = scoref_pair_tmp(i,3);
        label_pair{i} = sprintf('%s/%s',label_piano{p1},label_piano{p2});
    end
    
    score_avg_noise = median(score_pool_noise);
    scoreU_noise    = prctile(score_pool_noise,95);
    scoreL_noise    = prctile(score_pool_noise, 5);
        
    [label_pair idx] = sort(label_pair);
    score_avg_noise  = score_avg_noise(idx);
    scoreL_noise = scoreL_noise(idx);
    scoreU_noise = scoreU_noise(idx);
    
    [xx idx2]       = sort(score_avg_noise,'descend'); % second ordering criterion
    score_avg_noise = score_avg_noise(idx2);
    scoreL_noise = scoreL_noise(idx2);
    scoreU_noise = scoreU_noise(idx2);
    label_pair      = label_pair(idx2);
    
    errorU    = scoreU_noise - score_avg_noise;
    errorL    = score_avg_noise - scoreL_noise; 
    
    offsetx = 0.1;
    % figure;
    xdata = 1:length(score_avg_noise);
    errorbar(xdata+offsetx,score_avg_noise,errorL,errorU,'bs','LineWidth',2); grid on
    Xlabel('Piano pair',FontSize)
    Ylabel('Most similar pair, Score[\%]',FontSize)
    
    xlim([min(xdata)-0.5 max(xdata)+0.5]);
    ylim([0 103])
    set(gca,'XTick',xdata);
    h(end+1) = gcf;
    h(end) = Figure2paperfigureT(h(end),1,1,plotOpts);
    set(gca,'XTickLabel',label_pair);
    set(gca,'FontSize',FontSize);
    hold on;
    legend('silence','in noise, SNR = 0 dB')
    
    % name = 'triadic-three-participants-SNR';
    % Saveas(h(end),sprintf('%s%s',outputdir,name),'epsc')
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function yenv = il_get_envelope(insig,fs,fc)

if nargin < 3
    fc = 20;
end

yin = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv = filtfilt(b,a,yin);


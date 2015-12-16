function l0P392_20150331_my_plot_example
% function l0P392_20150331_my_plot_example
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirmain = 'D:\Documenten-TUe\09-Training+activities\Advanced-perception\2014-2015-Q3\Assignment-4\';
dirOut  = [dirmain 'Figures-new' delim]; 
dirData = [dirmain 'Assignment4_Instruction-Files' delim];

Mkdir(dirOut);

bDiary = 0;
Diary(mfilename,bDiary,dirOut);

bPart1 = 1;
bPart2 = 1;

if bPart1
    
[audio_delay, asynchronous, synchronous, N] = textread([dirData 'SJ2_1.csv'], '%d; %d; %d; %d', 'headerlines', 1);

% Plot options:
set(gcf,'DefaultLineLineWidth',3)
set(gcf,'DefaultAxesLineWidth',3)
set(gcf,'defaulttextfontsize',16)
set(gcf,'defaultaxesfontsize',16)
set(gca,'box','on')        
        
hold on;
plot(audio_delay, asynchronous./N, '-s', audio_delay,  synchronous./N, '-s');
xlabel('Audio delay (ms)');
ylabel('Response proportion');
axis([min(audio_delay)-25 max(audio_delay)+25 0.0 1.0]);
title('Raw SJ2 data for subject 1');
grid;
legend('Asynchronous', 'Synchronous', 'Location', 'East');
hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart2 % as sent by Armin

close all;
 
bPlotRaw = 1;
N_participants = 6;
hFigures = [];

%% First load the data sets:
% Filenames are SJ2_1.csv - SJ2_6.csv, SJ3_1.csv - SJ3_6.csv, TOJ_1.csv -
% TOJ_2.csv
for p = 1:6
    % Put data for different participants/data sets in different columns:
    [SJ2_audio_dly(:,p), SJ2_asynchro(:,p), SJ2_synchro(:,p), SJ2_N(:,p)]           = textread(sprintf('SJ2_%i.csv', p), '%d; %d; %d; %d', 'headerlines', 1);
    [SJ3_audio_dly(:,p), SJ3_audio_first(:,p), SJ3_synchro(:,p), SJ3_video_first(:,p), SJ3_N(:,p)]= textread(sprintf('SJ3_%i.csv', p), '%d; %d; %d; %d; %d', 'headerlines', 1);
    [TOJ_audio_dly(:,p), TOJ_audio_first(:,p), TOJ_video_first(:,p), TOJ_N(:,p)]    = textread(sprintf('TOJ_%i.csv', p), '%d; %d; %d; %d', 'headerlines', 1);
end 

% To simplify the understanding of this script:
N_trials = SJ2_N(1,1); % Be careful, I am assuming that SJ2_N, SJ3_N and TOJ_N are always 60
test_time_dly = SJ2_audio_dly(:,1);

%% 1.1 a) Plot the data
for p=1:N_participants
    
    if bPlotRaw
        % Plot SJ2:
        figure; % new figure for each participant & respone category
        % Plot options:
        Il_prepare_fig_default; % do this each time you create a figure
        
        subplot(3,1,1)
        
        hold on;
        plot(test_time_dly, SJ2_asynchro(:,p)./N_trials, '-s', SJ2_audio_dly(:,p), SJ2_synchro(:,p)./N_trials, '-s');
        xlabel('Audio delay (ms)');
        ylabel('Response proportion');
        axis([min(test_time_dly)-25 max(test_time_dly)+25 0.0 1.0]);
        title(sprintf('Raw SJ2 data for participant %i',p))
        grid;
        legend('Asynchronous', 'Synchronous', 'Location', 'East');
        hold off;

        % Plot SJ3:
        subplot(3,1,2)
        % Il_prepare_fig_default; % Uncomment if you use different figures

        hold on;
        plot(test_time_dly, SJ3_audio_first(:,p)./N_trials, '-s', SJ3_audio_dly(:,p), SJ3_synchro(:,p)./N_trials, '-s', SJ3_audio_dly(:,p), SJ3_video_first(:,p)./N_trials, '-s');
        xlabel('Audio delay (ms)');
        ylabel('Response proportion');
        axis([min(SJ3_audio_dly(:,p))-25 max(SJ3_audio_dly(:,p))+25 0.0 1.0]);
        title(sprintf('Raw SJ3 data for participant %i',p))
        grid;
        legend('Audio first', 'Synchronous', 'Video first', 'Location', 'East');
        hold off;

        % Plot TOJ:
        subplot(3,1,3)
        % Il_prepare_fig_default; % Uncomment if you use different figures   

        hold on;
        plot(test_time_dly, TOJ_audio_first(:,p)./N_trials, '-s', TOJ_audio_dly(:,p), TOJ_video_first(:,p)./N_trials, '-s');
        xlabel('Audio delay (ms)');
        ylabel('Response proportion');
        axis([min(test_time_dly)-25 max(test_time_dly)+25 0.0 1.0]);
        title(sprintf('Raw TOJ data for participant %i',p))
        grid;
        legend('Audio first', 'Video first', 'Location', 'East');
        hold off;
        
        Maximise_figure;
        hFigures(end+1) = gcf;
    end; % if(boPlotRaw)
    
end;
 
%% 1.1 b) Estimate PSE & 1.1 c) compute the synchrony window

% b-1) by taking the (mean of the) intersections (determined by linear interpolation)
% b-2) by taking the weighted average
% c-1) Synchrony window = the range of time delays for which the response 'synchronous' gives the highest score
% c-2) Use the synchrony window to get an additional estimate of the PSS
 
for p = 1:N_participants
    
    TOJ         = LinearIntersection(test_time_dly, TOJ_audio_first(:,p)./N_trials, TOJ_video_first(:,p)./N_trials); % b-1
    TOJ_PSE(p)  = TOJ(1); % b-1
    
    SJ2_wavg(p) = WAvg(test_time_dly, SJ2_synchro(:,p)); % b-2
    SJ3_wavg(p) = WAvg(test_time_dly, SJ3_synchro(:,p)); % b-2
    
    SJ2         = LinearIntersection(test_time_dly,    SJ2_asynchro(:,p)./N_trials,     SJ2_synchro(:,p)./N_trials);
    SJ3_left    = LinearIntersection(test_time_dly, SJ3_audio_first(:,p)./N_trials,     SJ3_synchro(:,p)./N_trials);
    SJ3_right   = LinearIntersection(test_time_dly, SJ3_video_first(:,p)./N_trials,     SJ3_synchro(:,p)./N_trials);
    
    SJ2_sync_window(p) = SJ2(2,1) - SJ2(1,1);        % c-1, it is a length value, difference between synchrony boundaries
    SJ3_sync_window(p) = SJ3_right(1) - SJ3_left(1); % c-1
    
    SJ2_PSS(p)  = (SJ2(1,1)+SJ2(2,1))/2;        % c-2
    SJ3_PSS(p)  = (SJ3_left(1)+SJ3_right(1))/2; % c-2
    
end

% Show results:
disp('SJ2_PSE:')
var2latex( SJ2_PSS )

disp('SJ3_PSE:')
var2latex( SJ3_PSS )

disp('TOJ_PSE:')
var2latex( TOJ_PSE )

disp('SJ2_sync_window:')
var2latex( SJ2_sync_window )

disp('SJ3_sync_window:')
var2latex( SJ3_sync_window )

disp('SJ2_wavg:')
var2latex( SJ2_wavg )

disp('SJ3_wavg:')
var2latex( SJ3_wavg )
 
%% 1.2 Comparison of data across subjects

% a) Compute the correlation between PSE estimates:
[r, p] = corrcoef(SJ3_PSS, SJ2_PSS); % r - correlation, p - p-value
r_SJ3_SJ2 = r(1,2);
p_SJ3_SJ2 = p(1,2);
[r, p] = corrcoef(SJ3_PSS, TOJ_PSE);
r_SJ3_TOJ = r(1,2);
p_SJ3_TOJ = p(1,2);
[r, p] = corrcoef(SJ2_PSS, TOJ_PSE);
r_SJ2_TOJ = r(1,2);
p_SJ2_TOJ = p(1,2);
 
disp('Correlation: r_SJ2_TOJ - r_SJ3_SJ2 - r_SJ3_TOJ')
var2latex([r_SJ2_TOJ r_SJ3_SJ2 r_SJ3_TOJ] )

% b) Compute pair-wise differences:
SJ2_SJ3 = SJ2_PSS - SJ3_PSS;
vSJ2_SJ3 = var(SJ2_SJ3); % variance

SJ2_TOJ = SJ2_PSS - TOJ_PSE;
vSJ2_TOJ = var(SJ2_TOJ); % variance

SJ3_TOJ = SJ3_PSS - TOJ_PSE;
vSJ3_TOJ = var(SJ3_TOJ); % variance

var2latex([SJ2_SJ3 vSJ2_SJ3; SJ2_TOJ vSJ2_TOJ; SJ3_TOJ vSJ3_TOJ])
% c) Formulate hyphoteses:
% %       CHECK ASSUMPTIONS
% [h1,p1] = vartest2(SJ2_PSE,SJ3_PSE,0.05,'both');    % h1 = 0 means: null hyphotesis can be accepted (it is likely)
%                                                     % h1 = 1 means: null hyphotesis can be rejected (returned p-value is less than 0.05)
%                                                     % null hyphotesis: variances are equal
% [h2,p2] = vartest2(SJ2_PSE,TOJ_PSE,0.05,'both');
% [h3,p3] = vartest2(SJ3_PSE,TOJ_PSE,0.05,'both');

[h11,p11] = vartest2(SJ2_SJ3,SJ2_TOJ,0.05,'both'); % different variance (expected)
[h12,p12] = vartest2(SJ2_SJ3,SJ3_TOJ,0.05,'both'); % different variance (expected)
[h13,p13] = vartest2(SJ2_TOJ,SJ3_TOJ,0.05,'both'); % same variance (large variance)

% d) Comparing the size of the synchrony windows:
[h21,p21] = ttest(SJ2_sync_window,SJ3_sync_window,0.05,'both'); % 0 = they have the same length

Save_all_figures(hFigures,dirOut);

disp('')

end

if bDiary
	diary off
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline functions
function Il_prepare_fig_default

    set(gcf,'DefaultLineLineWidth',3)
    set(gcf,'DefaultAxesLineWidth',3)
    set(gcf,'defaulttextfontsize',16)
    set(gcf,'defaultaxesfontsize',16)
    set(gca,'box','on')        

end
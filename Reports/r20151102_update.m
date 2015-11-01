function r20151102_update
% function r20151102_update
%
% 1. Description:
%       Plotting results obtained from r20151016_update_dau1997_jepsen2008.m
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 31/10/2015
% Last update on: 31/10/2015 
% Last use on   : 31/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

dir_out = [Get_TUe_paths('outputs') '20151102-validation' delim];
Mkdir(dir_out);

% CASP = model 103 (sigma = 0.59)
% PEMO = model 101 (sigma = 0.68; modfiltertype = 'dau1997wLP')

% TODO:
%       Intensity discrimination: run with other values of sigma
%       TMTF: re-run jepsen2008, 31-Hz NBN

close all
hFig4 = []; % handles for the generated figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Results of intensity discrimination task (similar to Fig 3 in Jepsen2008)
%%% Audio files:
% 1.1 Intensity discrimination with tones:
%   'AMTControl-examples/tone-f-1000-Hz-at-60-dB-dur-800-ms.wav'
%   'AMTControl-examples/tone-f-1000-Hz-at-42-dB-dur-800-ms.wav'
% 1.2 Intensity discrimination with BBN:
%   'audio-20151006/jepsen2008-BW-at-60-dB-dur-500-ms-2.wav'
%   'audio-20151006/jepsen2008-BW-at-42-dB-dur-500-ms-2.wav'
    
% raw1 - intensity discrimination with 1-kHz tones
% raw2 - intensity discrimination with BBN (100 Hz - 8 kHz)
% raw3_j2008 - intensity discrimination with 1-kHz tones, different sigmas (only jepsen2008)
% raw4_j2008 - intensity discrimination with BBN, different sigmas (only jepsen2008)

testlevels = [20 30 40 50 60 70];
minL = 20; maxL = 70;
minY = 0.3; maxY = 1.5;
     
raw1_d1997 = [1.2778    1.0299    0.8279    0.7416    0.6639    0.5941]; % sigma = 0.68
raw2_d1997 = [0.5941    0.5314    0.6639    0.7623    0.8279    0.8279]; % sigma = 0.68

raw1_j2008 = [0.8279    0.6639    0.5314    0.5314    0.5314    0.5314]; % sigma = 0.59
raw2_j2008 = [0.6639    0.5314    0.5941    0.6639    0.8279    1.0299]; % sigma = 0.59
               
figure;
subplot(1,2,1)
plot(testlevels, raw1_d1997,'rs-','LineWidth',2); hold on
plot(testlevels, raw1_j2008,'ko-','LineWidth',2); grid on
xlabel('Standard level [dB]'); xlim([minL maxL])
ylabel('\Delta at threshold [dB]'); ylim([minY maxY])
title('A. test 1-kHz tones');
set(gca,'XTick',testlevels);
set(gca,'YTick',[minY:0.1:maxY])
     
subplot(1,2,2)
plot( testlevels, raw2_d1997,'rs-','LineWidth',2); hold on
plot( testlevels, raw2_j2008,'ko-','LineWidth',2); grid on
xlabel('Standard level [dB]'); xlim([minL maxL])
ylabel('\Delta at threshold [dB]'); ylim([minY maxY])
title('B. test BBNs');
set(gca,'XTick',testlevels);
set(gca,'YTick',[minY:0.1:maxY])
     
legend('PEMO', 'CASP','Location','NorthWest');

hFig4(end+1) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The results of the CASP model with different sigma values:
% raw3 - results for sine tones
% raw4 - results for BBNs

sigmas      = [0.5 0.55 0.59 0.65 0.7 0.8];
             %[0.5 0.6 0.615 0.63 0.7 0.85]; % sigmas in old version

LineColor  = {'bs-';'m>-.';'ko-';'r--';'g>-';'r<-'};
LineWidth  = [1.5 1 2 1 1 2];


raw3_j2008 = [0.6639    0.5314    0.4752    0.4248    0.4752    0.4752; ... % sigma = 0.5
              0.7416    0.5941    0.5314    0.4752    0.5314    0.4752; ... % sigma = 0.55
              0.8279    0.6639    0.5314    0.5314    0.5314    0.5314; ... % sigma = 0.59
              0.9237    0.6639    0.5941    0.5941    0.5941    0.5941; ... % sigma = 0.65
              0.9237    0.7416    0.6639    0.5941    0.6639    0.6639; ... % sigma = 0.7
              1.1476    0.8279    0.7416    0.7416    0.7416    0.7416];    % sigma = 0.8

raw4_j2008 = [0.5941    0.4248    0.5314    0.5314    0.6639    0.8279; ... % sigma = 0.5
              0.5941    0.4752    0.5314    0.5941    0.7416    0.9237; ... % sigma = 0.55
              0.6639    0.5314    0.5941    0.6639    0.8279    1.0299; ... % sigma = 0.59
              0.7416    0.5941    0.6639    0.7416    0.9237    1.1476; ... % sigma = 0.65
              0.8279    0.5941    0.7416    0.7416    0.9237    1.1476; ... % sigma = 0.7
              0.9237    0.7416    0.8279    0.8279    1.0299    1.4216];    % sigma = 0.8

figure;
for i = 1:length(LineColor)
    subplot(1,2,1)
    plot(testlevels, raw3_j2008(i,:), LineColor{i},'LineWidth',LineWidth(i)), hold on, grid on
    if i==1;
        title('C. test 1-kHz tones');
        xlabel('Standard level [dB]')
        ylabel('\Delta at threshold [dB]')
        xlim([minL maxL])
        ylim([minY maxY])
        set(gca,'XTick',testlevels);
        set(gca,'YTick',[minY:0.1:maxY]);
    end
    subplot(1,2,2)
    plot(testlevels, raw4_j2008(i,:), LineColor{i},'LineWidth',LineWidth(i)), hold on, grid on
    if i==1;
        title('D. test BBNs');
        xlabel('Standard level [dB]')
        ylabel('\Delta at threshold [dB]')
        xlim([minL maxL])
        ylim([minY maxY])
        set(gca,'XTick',testlevels);
        set(gca,'YTick',[minY:0.1:maxY]);
    end
    Leg{i} = ['\sigma = ' num2str(sigmas(i))];
end
legend(Leg,'Location','NorthWest')

hFig4(end+1) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Temporal modulation transfer functions (similar to Fig 8 in Jepsen2008)

% Results correspond to 3-Hz and 31-Hz data as obtained on 30/10/2015 for Dau1997
% and 3-Hz and 31-Hz data of the second run of Jepsen2008

% raw1 - TMTF,  3-Hz NBN
% raw2 - TMTF, 31-Hz NBN

minY = -35;
maxY = 0;
fmodtest = [3 5 10 20 50 100 120];

raw1_d1997 = [-11.5000  -20.2500  -30.0000  -28.5000  -27.5000  -27.7500  -26.2500]; % supra level = -6 dB
%%% 3-Hz, Dau1997. supraT level = -6 dB (run on 01/11/2015)
r1_d75 = [-10.7500  -19.0000  -28.0000  -27.0000  -26.0000  -26.5000  -25.0000]; r1_d75 = r1_d75 - raw1_d1997;
r1_d25 = [-13.2500  -21.7500  -31.7500  -29.2500  -30.0000  -28.7500  -29.0000]; r1_d25 = raw1_d1997 - r1_d25;

raw2_d1997 = [-5.5000   -9.5000  -13.7500  -12.7500  -19.2500  -18.0000  -19.2500]; % supra level = -1 dB
%%% 31-Hz, Dau1997. supraT level = -1 dB (run on 01/11/2015)
r2_d75 = [-4.5000   -7.5000  -12.2500  -12.0000  -18.0000  -18.0000  -16.7500]; r2_d75 = r2_d75 - raw2_d1997;
r2_d25 = [-6.0000  -10.0000  -14.0000  -14.0000  -20.5000  -19.0000  -22.7500]; r2_d25 = raw2_d1997 - r2_d25;

raw1_j2008 = [-11.0000  -13.0000  -22.7500  -22.7500  -18.0000  -22.0000  -16.7500]; 
%%% 3-Hz, Jepsen2008. supraT level = -6 dB (run on 01/11/2015)
r1_j75 = [ -8.2500  -11.0000  -21.2500  -22.0000  -17.5000  -20.0000  -14.5000]; r1_j75 = r1_j75 - raw1_j2008;
r1_j25 = [-13.0000  -15.2500  -24.2500  -26.5000  -19.5000  -24.7500  -18.2500]; r1_j25 = raw1_j2008 - r1_j25;

raw2_j2008 = [-1.7500   -6.5000  -11.2500  -11.2500  -17.5000  -17.2500  -17.0000]; 
%%% 31-Hz, Jepsen2008. supraT level = -1 dB, first fmod might be biased by 0 response (run on 01/11/2015)
r2_j75 = [-1.5000   -5.7500  -10.2500   -9.2500  -16.2500  -16.2500  -16.0000]; r2_j75 = r2_j75 - raw2_j2008;
r2_j25 = [-2.5000   -7.5000  -11.7500  -14.2500  -18.2500  -18.5000  -17.0000]; r2_j25 = raw2_j2008 - r2_j25;

figure;
subplot(1,2,1)
errorbar(1:length(fmodtest), raw1_d1997,r1_d25, r1_d75, 'rs-','LineWidth',2); hold on
errorbar(1:length(fmodtest), raw1_j2008,r1_j25, r1_j75,'ko-','LineWidth',2); grid on
xlabel('Modulation frequency [Hz]'); xlim([0 length(fmodtest)+1])
ylabel('Modulation depth [dB]'); ylim([minY maxY])
title('A. fmod = 3 Hz');
set(gca,'XTick',1:length(fmodtest));
set(gca,'XTickLabel',fmodtest);
set(gca,'YTick',[minY:5:maxY])
legend('PEMO', 'CASP','Location','NorthWest');

subplot(1,2,2)
errorbar(1:length(fmodtest), raw2_d1997,r2_d25,r2_d75,'rs-','LineWidth',2); hold on
errorbar(1:length(fmodtest), raw2_j2008,r2_j25,r2_j75,'ko-','LineWidth',2); grid on
xlabel('Modulation frequency [Hz]'); xlim([0 length(fmodtest)+1])
ylabel('Modulation depth [dB]'); ylim([minY maxY])
title('B. fmod = 31 Hz');
set(gca,'XTick',1:length(fmodtest));
set(gca,'XTickLabel',fmodtest);
set(gca,'YTick',[minY:5:maxY])
     
hFig4(end+1) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Forward masking: on- and off-frequency listening (similar to Fig 7 in Jepsen2008)

testlevels = [  40 60 70 80; ...
                60 70 80 85];
minL = min(min(testlevels));
maxL = max(max(testlevels));
     
raw1_d1997 = [26.2500   39.0000   50.5000   63.5000; ... % on-freq,  0-ms
              20        34        43        52     ];    % on-freq, 30-ms
r1_d75     = [26.5000   39.0000   50.5000   64.0000; ... % on-freq,  0-ms
              20        34        43        52     ];    % on-freq, 30-ms
r1_d25     = [26.0000   39.0000   50.2500   63.0000; ... % on-freq,  0-ms
              20        34        43        52     ];    % on-freq, 30-ms

raw2_d1997 = [14.0000   26.0000   37.0000   42.0000; ... % off-freq,  0-ms
               8        12        17        21    ];     % off-freq, 30-ms
r2_d75     = [14.0000   26.0000   37.0000   42.0000; ... % off-freq,  0-ms
               8        12        17        21    ];     % off-freq, 30-ms              
r2_d25     = [14.0000   25.5000   37.0000   42.0000; ... % off-freq,  0-ms
               8        12        17        21    ];     % off-freq, 30-ms

raw1_j2008 = [27.0000   30.0000   33.0000   38.0000; ... % on-freq,  0-ms
              25        30        32        35    ];     % on-freq, 30-ms
r1_j75     = [27.0000   30.0000   33.0000   38.0000; ... % on-freq,  0-ms
              25        30        32        35    ];     % on-freq, 30-ms
r1_j25     = [27.0000   30.0000   33.0000   38.0000; ... % on-freq,  0-ms
              25        30        32        35    ];     % on-freq, 30-ms

raw2_j2008 = [21.0000   27.0000   34.2500   43.0000; ... % off-freq,  0-ms
              18        20        25        29    ];     % off-freq, 30-ms
r2_j75     = [21.0000   27.0000   34.7500   43.0000; ... % off-freq,  0-ms
              18        20        25        29    ];     % off-freq, 30-ms
r2_j25     = [21        27        34        43     ; ... % off-freq,  0-ms
              18        20        25        29    ];     % off-freq, 30-ms

figure;
subplot(1,2,1)
plot(testlevels(1,:), raw1_d1997(1,:),'r>-','LineWidth',2); hold on
plot(testlevels(1,:), raw1_d1997(2,:),'ro--');
plot(testlevels(1,:), raw1_j2008(1,:),'k>-' );
plot(testlevels(1,:), raw1_j2008(2,:),'ko--'); grid on
xlabel('Masker level [dB SPL]'); xlim([minL maxL])
ylabel('Threshold [dB SPL]')
title('A. On-freq.');

subplot(1,2,2)
plot( testlevels(2,:), raw2_d1997(1,:),'r>-','LineWidth',2); hold on
plot( testlevels(2,:), raw2_d1997(2,:),'ro--')
plot( testlevels(2,:), raw2_j2008(1,:),'k>-' )
plot( testlevels(2,:), raw2_j2008(2,:),'ko--'); grid on
xlabel('Masker level [dB SPL]'); xlim([minL maxL])
ylabel('Threshold [dB SPL]')
title('B. Off-freq.');

legend('PEMO,  0-ms', ...
       'PEMO, 30-ms', ...
       'CASP,  0-ms', ...
       'CASP, 30-ms','Location','NorthWest');

hFig4(end+1) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Forward masking: tone-in-tone

% jepsen2008:
    %   -var = 0.59, gain_after_drnl = 17.86 dB
    
toffsetonset = [-20 -10 -5 0 10 20 30]; % corresponding to ([480 490 495 500:10:530])*1e-3;
testlevels = [40; 60; 80];

TTh_j2008 = [24.0000   29.0000   28.4375   27.5000   27.5000   26.2500   25.0000; ...
             55.4375   57.0000   36.0000   30.0000   30.0000   31.2500   30.0000; ...
             76.1875   78.1250   54.9375   37.5000   36.2500   36.2500   35.0000];
Prct75_j2008 = [26.2500   30.1250   30.0000   27.5000   27.5000   26.2500   25.0000; ...
                56.9375   57.6875   36.2500   30.0000   30.0000   31.2500   30.0000; ...
                78.4375   79.3750   55.0000   37.5000   36.2500   36.2500   35.0000];
Prct25_j2008 = [23.3125   25.8750   27.8125   27.5000   27.5000   26.2500   25.0000; ...
                54.4375   56.2500   35.6250   30.0000   30.0000   31.2500   30.0000; ...
                74.9375   74.5625   54.3125   37.5000   36.2500   36.2500   35.0000];

errorL_j2008 = TTh_j2008-Prct25_j2008;
errorU_j2008 = Prct75_j2008-TTh_j2008;

TTh_d1997 =   [ 25.5000   32.0000   31.2500   26.2500   25.0000   22.5000   20.0000; ...
                49.3750   53.9375   50.7500   39.7500   40.0000   37.5000   33.7500; ...
                76.3125   71.5625   71.8750   63.7500   58.7500   56.2500   51.2500];
Prct75_d1997 = [29.5000   34.0625   31.5625   26.2500   25.0000   22.5000   20.0000; ...
                53.1250   55.1875   51.3750   40.0000   40.0000   37.5000   33.7500; ...
                77.2500   74.6875   72.5000   63.7500   58.7500   56.2500   51.2500];
Prct25_d1997 = [24.0000   30.4375   31.1250   26.2500   25.0000   22.5000   20.0000; ...
                48.7500   53.2500   49.2500   39.3750   40.0000   37.5000   33.7500; ...
                73.6875   70.6250   70.8750   63.7500   58.7500   56.2500   51.2500];

errorL_d1997 = TTh_d1997-Prct25_d1997;
errorU_d1997 = Prct75_d1997-TTh_d1997;

figure;
subplot(1,2,1)
errorbar(toffsetonset,TTh_j2008(3,:),errorL_j2008(3,:),errorU_j2008(3,:),'ko-'); grid on; hold on
errorbar(toffsetonset,TTh_j2008(2,:),errorL_j2008(2,:),errorU_j2008(2,:),'k>-');
errorbar(toffsetonset,TTh_j2008(1,:),errorL_j2008(1,:),errorU_j2008(1,:),'ks-'); 
xlabel('Offset-onset interval [ms]');
ylabel('Masked threshold [dB SPL]')
legend(sprintf('M at %.0f dB',testlevels(3)), ...
       sprintf('M at %.0f dB',testlevels(2)), ...
       sprintf('M at %.0f dB',testlevels(1)) );
ha = gca;
title('Forward masking (jepsen2008)');

subplot(1,2,2)
errorbar(toffsetonset,TTh_d1997(3,:),errorL_d1997(3,:),errorU_d1997(3,:),'ro-'); grid on; hold on
errorbar(toffsetonset,TTh_d1997(2,:),errorL_d1997(2,:),errorU_d1997(2,:),'r>-');
errorbar(toffsetonset,TTh_d1997(1,:),errorL_d1997(1,:),errorU_d1997(1,:),'rs-'); 
xlabel('Offset-onset interval [ms]');
ylabel('Masked threshold [dB SPL]')
legend(sprintf('M at %.0f dB',testlevels(3)), ...
       sprintf('M at %.0f dB',testlevels(2)), ...
       sprintf('M at %.0f dB',testlevels(1)) );
ha(end+1) = gca;
title('Forward masking (dau1997)');

linkaxes(ha,'xy');
    
hFig4(end+1) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Forward masking: off-frequency listening

toffsetonset = [-20 -10 -5 0 10 20 30]; % corresponding to ([480 490 495 500:10:530])*1e-3;
testlevels = [40; 60; 80];

% Jepsen2008, off frequency
TTh_j2008 = [13.7500   13.7500   15.0000   15.0000   15.0000   15.0000   15.0000; ...
             24.7500   24.8750   22.5000   21.2500   18.7500   17.5000   17.5000; ...
             53.7500   52.6250   42.5000   34.2500   28.7500   27.5000   25.0000];
Prct75_j2008 = [13.7500   13.7500   15.0000   15.0000   15.0000   15.0000   15.0000; ...
                25.0000   25.0000   22.5000   21.2500   18.7500   17.5000   17.5000; ...
                54.0000   53.1875   42.5000   34.2500   28.7500   27.5000   25.0000];
Prct25_j2008 = [13.7500   13.7500   15.0000   15.0000   15.0000   15.0000   15.0000; ...
                24.3750   24.6250   22.5000   21.2500   18.7500   17.5000   17.5000; ...
                53.7500   52.5625   42.5000   34.2500   28.7500   27.5000   25.0000];

errorL_j2008 = TTh_j2008-Prct25_j2008;
errorU_j2008 = Prct75_j2008-TTh_j2008;

TTh_d1997 =   [  5.0000    5.0000    5.0000    5.0000    5.0000    5.0000    5.0000; ...
                11.2500   11.2500   17.1250   13.7500   10.0000    8.7500    7.5000; ...
                31.2500   31.2500   37.5625   36.2500   25.0000   21.2500   17.5000];
Prct75_d1997 = [ 5.0000    5.0000    5.0000    5.0000    5.0000    5.0000    5.0000; ...
                11.2500   11.2500   17.3750   13.7500   10.0000    8.7500    7.5000; ...
                31.2500   31.2500   37.6250   36.2500   25.0000   21.2500   17.5000];
Prct25_d1997 = [ 5.0000    5.0000    5.0000    5.0000    5.0000    5.0000    5.0000; ...
                11.2500   11.2500   16.6875   13.7500   10.0000    8.7500    7.5000; ...
                31.2500   31.2500   37.5000   36.2500   25.0000   21.2500   17.5000];

errorL_d1997 = TTh_d1997-Prct25_d1997;
errorU_d1997 = Prct75_d1997-TTh_d1997;

figure;
subplot(1,2,1)
errorbar(toffsetonset,TTh_j2008(3,:),errorL_j2008(3,:),errorU_j2008(3,:),'ko-'); grid on; hold on
errorbar(toffsetonset,TTh_j2008(2,:),errorL_j2008(2,:),errorU_j2008(2,:),'k>-');
errorbar(toffsetonset,TTh_j2008(1,:),errorL_j2008(1,:),errorU_j2008(1,:),'ks-'); 
xlabel('Offset-onset interval [ms]');
ylabel('Masked threshold [dB SPL]')
legend(sprintf('M at %.0f dB',testlevels(3)), ...
       sprintf('M at %.0f dB',testlevels(2)), ...
       sprintf('M at %.0f dB',testlevels(1)) );
ha = gca;
title('Forward masking, off-frequency masker (jepsen2008)');

subplot(1,2,2)
errorbar(toffsetonset,TTh_d1997(3,:),errorL_d1997(3,:),errorU_d1997(3,:),'ro-'); grid on; hold on
errorbar(toffsetonset,TTh_d1997(2,:),errorL_d1997(2,:),errorU_d1997(2,:),'r>-');
errorbar(toffsetonset,TTh_d1997(1,:),errorL_d1997(1,:),errorU_d1997(1,:),'rs-'); 
xlabel('Offset-onset interval [ms]');
ylabel('Masked threshold [dB SPL]')
legend(sprintf('M at %.0f dB',testlevels(3)), ...
       sprintf('M at %.0f dB',testlevels(2)), ...
       sprintf('M at %.0f dB',testlevels(1)) );
ha(end+1) = gca;
title('Forward masking, off-frequency masker (dau1997)');

linkaxes(ha,'xy');
ylim([2 60]);

hFig4(end+1) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Signal integration:
testdurs = [10 20 40 70 150]; % in ms
% Jepsen2008: (25 avgs., 10-ms ramps)
raw_j2008 = [66.2500   60.5000   58.5000   58.5000   55.7500];
Std75 = [67.5000   60.7500   58.7500   59.5000   56.7500];
Std25 = [65.2500   60.2500   58.2500   57.2500   55.5000];

rj75 = Std75 - raw_j2008;
rj25 = raw_j2008 - Std25;

% Dau1997 % (25 avgs., 10-ms ramps)
raw_d1997 = [61.2500   58.7500   59.0000   57.2500   58.5000];
Std75 = [61.5000   60.2500   59.0000   57.7500   59.0000];
Std25 = [61.0000   58.0000   58.0000   56.5000   57.7500];
rd75 = Std75 - raw_d1997;
rd25 = raw_d1997 - Std25;

minY = 50;
maxY = 75;

xoffset = 0.05;
xtestdurs = 1:length(testdurs);
figure;
errorbar(xtestdurs-xoffset, raw_d1997,rd25,rd75,'rs-','LineWidth',2); hold on
errorbar(xtestdurs+xoffset, raw_j2008,rj25,rj75,'ko-','LineWidth',2); grid on
xlabel('Tone duration [ms]'); xlim([0 max(xtestdurs)+1])
ylabel('Masked threshold [dB SPL]'); ylim([minY maxY])
title('Signal integration experiment');

set(gca,'XTick',xtestdurs);
set(gca,'XTickLabel',testdurs);
set(gca,'YTick',[minY:5:maxY])
     
legend('PEMO', 'CASP','Location','NorthWest');

hFig4(end+1) = gcf;

Save_all_figures(hFig4,dir_out);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

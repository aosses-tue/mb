function r20150801_understanding_IQR
% function r20150801_understanding_IQR
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 02/08/2015
% Last update on: 02/08/2015 % Update this date manually
% Last use on   : 02/08/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bDiary = 0;
Diary(mfilename,bDiary);

bDoSRT_Leuven = 0;
bDoPooledPilot = 1;

if bDoSRT_Leuven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   25% of 24 =  6 responses
%   75% of 24 = 18 responses

% Strat = 10
% Subject, strat; SRT 
data = [11	10	-2.33; ...
 11	10	-5   ; ...
 11	10	-2.33; ...
13	10	-1.67; ...
13	10	-4.33; ...
14	10	-1   ; ...
14	10	-1	 ; ...
15	10	-3	 ; ...
15	10	-4.33; ...
16	10	 2.33; ...
16	10	1; ...
16	10	1; ...
17	10	1; ...
17	10	-0.33; ...
17	10	3];
SRT1 = data(:,3);

% Strat = 0
data = [11 	0	-1.67; ...
11	0	-3.67; ...
12 	0	1; ...
12	0	-1.67; ...
12	0	2.33; ...
13 	0	-1.67; ...
13	0	-2.33; ...
14	0	-1; ...
14	0	2.33; ...
15	0	-2.33; ...
15	0	-1; ...
16	0	-1; ...
16	0	2.33; ...
16	0	5; ...
16	0	0.33
17	0	2.33; ...
17	0	1];
SRT2 = data(:,3);

% Strat = 1
data = [11	1	-2.33; ...
11	1	-3.67; ...
12	1	3; ...
12	1	5.67; ...
12	1	3.67; ...
13	1	-3; ...
13	1	-1; ...
14	1	-0.33; ...
14	1	1; ...
15	1	-3; ...
15	1	-3; ...
16	1	3.67; ...
16	1	2.33; ...
16	1	7; ...
16	1	-1; ...
17	1	1.67; ...
17	1	3.67; ...
17	1	3.67; ...
17	1	3; ...
17	1	0.33; ...
17	1	2.33];
SRT3 = data(:,3);

perc = [percentile(SRT1,50) percentile(SRT1,25) percentile(SRT1,75); ...
        percentile(SRT2,50) percentile(SRT2,25) percentile(SRT2,75); ...
        percentile(SRT3,50) percentile(SRT3,25) percentile(SRT3,75)];

meanstd = [ mean(SRT1) std(SRT1); ...
            mean(SRT2) std(SRT2); ...
            mean(SRT3) std(SRT3)];

[Mean(1,1),errorL(1,1),errorU(1,1)] = Prepare_errorbar_perc(SRT1,25,75);
[Mean(2,1),errorL(2,1),errorU(2,1)] = Prepare_errorbar_perc(SRT2,25,75);
[Mean(3,1),errorL(3,1),errorU(3,1)] = Prepare_errorbar_perc(SRT3,25,75);
              
figure; 
errorbar([1 2 3]-0.10,Mean,errorL,errorU,'bx'); hold on
errorbar([1 2 3]+0.10,meanstd(:,1),meanstd(:,2),'ro');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bDoPooledPilot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirres = '/home/alejandro/Documenten/Documenten-TUe/09-Training+activities/Master-thesis/01-2015-Rodrigo/20150729-pilot-APEX/';
dirres = '/home/alejandro/Documenten/Documenten-TUe/02-Experiments/2015-APEX-Rodrigo/APEX_shared/experiment/fluctuation_strength_results/';
fileres = { 'AM-fm-S01-RO.apr', ...
            'AM-fm-S02-CH.apr', ...
            'AM-fm-S03-ED.apr', ...
            'AM-fm-S04-AN.apr', ...
            'AM-fm-S05-BE.apr', ...
            'AM-fm-S06-RA.apr', ...
            'AM_tones-fm-results-AO-test.apr', ...
            'AM_tones-fm-slider-results-AO-test.apr'};

% fileres = { 'AM-fm-S01-RO.apr', ...
%             'AM-fm-S05-BE.apr', ...
%             'AM-fm-S06-RA.apr', ...
%             'AM_tones-fm-results-AO-test.apr', ...
%             'AM_tones-fm-slider-results-AO-test.apr'};
   
Nfiles = length(fileres);
dataraw1 = nan(4*Nfiles,11);
dataraw2 = nan(4*Nfiles,11);

for i = 1:length(fileres)
    
    resTmp = r20150802_plot_FS_results_v2([dirres fileres{i}]);
    Ntrials = size(resTmp.scores1norm,2); % 'local' number of trials: 9 or 11
    dataraw1(4*(i-1)+1:4*i,1:Ntrials) = resTmp.scores1norm;
    dataraw2(4*(i-1)+1:4*i,1:Ntrials) = resTmp.scores2norm;
    
end

[m1 s1 ci1 iq1] = Get_mean(dataraw1);
[m2 s2 ci2 iq2] = Get_mean(dataraw2);

test_fmod = [0 0.25 0.5 1 2 4 8 16 32 64 128];
text_XLabel = 'f_m_o_d [Hz]';
text_YLabel = 'Relative FS [%]';
Ntest = length(test_fmod);

figure; 
subplot(2,1,1)
errorbar([1:Ntest]-0.10,iq1(1,:),iq1(2,:),iq1(3,:),'bx'); hold on
errorbar([1:Ntest]+0.10,m1,s1,'ro');
xlabel(text_XLabel)
ylabel(text_YLabel)
grid on;
xlim([0.5 Ntest+0.5])
ha = gca;
set(ha,'XTickLabel',test_fmod)
title('FS using Standard1')
legend('IQR','std')

subplot(2,1,2)
errorbar([1:Ntest]-0.10,iq2(1,:),iq2(2,:),iq2(3,:),'bx'); hold on
errorbar([1:Ntest]+0.10,m2,s2,'ro');
xlabel(text_XLabel)
ylabel(text_YLabel)
grid on;
xlim([0.5 Ntest+0.5])
ha = gca;
set(ha,'XTickLabel',test_fmod)
title('FS using Standard2')
legend('IQR','std')

%%%

figure; 
errorbar([1:Ntest]-0.10,iq1(1,:),iq1(2,:),iq1(3,:),'bx'); hold on
errorbar([1:Ntest]+0.10,iq2(1,:),iq2(2,:),iq2(3,:),'ro');
plot([1:Ntest],mean([iq1(1,:); iq2(1,:)]),'k>--','LineWidth',2)
xlabel(text_XLabel)
ylabel(text_YLabel)
grid on;
xlim([0.5 Ntest+0.5])
ha = gca;
set(ha,'XTickLabel',test_fmod)
title('Measured FS')
legend('std1','std2','combined')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
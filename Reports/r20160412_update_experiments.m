function r20160412_update_experiments
% function r20160412_update_experiments
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/04/2016
% Last update on: 12/04/2016 
% Last use on   : 12/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
dir_where = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\'; %Stage-7-Pilot-HW'

h = [];
FontSize = 14;
offset = 0.05;

dir_1 = [dir_where 'Stage-6-Pilot-FD+AO\Exp-2-constant\'];
files = Get_filenames(dir_1,'*FD*apr');
%     (1,1): piano_multi-constant-FD.apr
[perc proc] = quick_constant_multi([dir_1 files{1,1}]);

disp('')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_HW = [dir_where 'Stage-7-Pilot-HW' delim]; % Huihui
SNR = [-5 5 20];
Scores = nan(4,length(SNR));

files = Get_filenames(dir_HW,'*20-dB*apr');
% Expected output:
%     (1,1): piano_F3-GH-JBS-constant-05-HW-SNR-p20-dB.apr
%     (1,2): piano_multi-F3-C4-same-Piano-constant-HW-SNR-p20-dB.apr

[perc proc] = quick_constant_multi([dir_HW files{1,2}]);
idx = [3 4 1 2]; % 4t2t3, 4t3t2, 5t2t4, 5t4t2
perc = perc(idx);

Scores(:,3) = transpose(perc);

% double	1 x 8	procedure5_t2_t4
% double	1 x 8	procedure5_t4_t2
% double	1 x 8	procedure4_t2_t3
% double	1 x 8	procedure4_t3_t2

files = Get_filenames(dir_HW,'*p05-dB*apr');
[perc proc] = quick_constant_multi([dir_HW files{1,1}]);
idx = [1 3 2 4]; % 4t2t3, 4t3t2, 5t2t4, 5t4t2
perc = perc(idx);

Scores(:,2) = transpose(perc);

figure;
plot((1:length(SNR))-offset,Scores(1,:),'bs--'); hold on
plot((1:length(SNR))+offset,Scores(2,:),'bo--');
plot((1:length(SNR))-offset,Scores(3,:),'rs-','LineWidth',2);
plot((1:length(SNR))+offset,Scores(4,:),'ro-','LineWidth',2);

xlim([1 3.33])
ylim([0 103])
grid on
ha = gca;
xres = [2-0.67 2 length(SNR)];
set(ha,'XTick',xres);
set(ha,'XTickLabel',SNR);
set(ha,'FontSize',FontSize);
Xlabel('SNR [dB]',FontSize);
Ylabel('Scores [\%]',FontSize);

legend({'C_4: JBS51, t2 (test) / t3 (ref)', ...
        'C_4: JBS51, t3 (test) / t2 (ref)', ...
        'F_3: JBS50, t2 (test) / t4 (ref)', ...
        'F_3: JBS50, t4 (test) / t2 (ref)'},'FontSize',10,'Location','NorthWest')

h(end+1) = gcf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR = [-5 5 20];
ScoresInter = nan(1,3);
files = Get_filenames(dir_HW,'*GH*20-dB*apr');
% Expected output:
%     (1,1): piano_F3-GH-JBS-constant-05-HW-SNR-p20-dB.apr
%     (1,2): piano_multi-F3-C4-same-Piano-constant-HW-SNR-p20-dB.apr

[perc proc] = quick_constant_multi([dir_HW files{1,1}]);

ScoresInter(:,3) = transpose(perc);

% double	1 x 8	procedure5_t2_t4
% double	1 x 8	procedure5_t4_t2
% double	1 x 8	procedure4_t2_t3
% double	1 x 8	procedure4_t3_t2

files = Get_filenames(dir_HW,'*GH*05-dB*apr');
[perc proc] = quick_constant_multi([dir_HW files{1,1}]);

ScoresInter(:,1) = transpose(perc);

% figure;
plot(xres,ScoresInter,'k>','LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files = Get_filenames(dir_HW,'*adaptive*apr');
[o1 o2 o3 o4] = quick_staircases([dir_HW files{1,1}]);

N4ave = 8;

figure;
Get_mAFC_reversals(o4.staircases.procedure05);
reve(1,:) = Get_mAFC_reversals(o4.staircases.procedure05);
Xlabel('Trial number',FontSize)
Ylabel('Adaptive parameter, SNR [dB]',FontSize)
h(end+1) = gcf;
grid on
title(sprintf('Proc 05: Threshold = %.2f [dB]',median(reve(1,end-N4ave+1:end)))) 


figure;
Get_mAFC_reversals(o4.staircases.procedure50);
reve(2,:) = Get_mAFC_reversals(o4.staircases.procedure50);
Xlabel('Trial number',FontSize)
Ylabel('Adaptive parameter, SNR [dB]',FontSize)
h(end+1) = gcf;
grid on
title(sprintf('Proc 50: Threshold = %.2f [dB]',median(reve(2,end-N4ave+1:end))))

for i = 1:length(h)
    Saveas(h(i),sprintf('results-HW-%.0f',i),'epsc');
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

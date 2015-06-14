function r20150612_update
% function r20150612_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 12/06/2015
% Last update on: 12/06/2015 % Update this date manually
% Last use on   : 12/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

close all

f       = [125  250  500 1000 2000 4000 8000];
lvl100  = [NaN  NaN  NaN 104.9 NaN  NaN 100.9; ...
           NaN  NaN  NaN 104.6 NaN  NaN 102.3];
lvl050  = [NaN  NaN  NaN  99.0 NaN  NaN  95.0 ; ...
           NaN  NaN  NaN  98.6 NaN  NaN  96.4];
lvl010  = [88.6 89.1 85.5 85.1 87.2 88.2 77.6; ...
           88.6 89.0 84.5 85.0 86.9 88.2 81.6]; % Avg. 2 measurements

       
Table1 = [f f; lvl100(1,:) lvl100(2,:); lvl050(1,:) lvl050(2,:); lvl010(1,:) lvl010(2,:)];
var2latex(Table1);
       
vtg100  = [0.707	0.725	0.739	0.743	0.743	0.699	0.542; ...
           0.705	0.722	0.737	0.742	0.74	0.698	0.541] * 1000;
vtg050  = [0.359	0.365	0.372	0.374	0.375	0.376	0.273; ...
           0.356	0.365	0.374	0.375	0.376	0.354	0.274] * 1000;
vtg010  = [0.0716	0.073	0.0743	0.0746	0.0748	0.075	0.076; ...
           0.0714	0.0734	0.075	0.076	0.0754	0.0756	0.0767] *1000;

Table2 = [f f; vtg100(1,:) vtg100(2,:); vtg050(1,:) vtg050(2,:); vtg010(1,:) vtg010(2,:)];
var2latex(Table2);

Sens = 104.75 - 20*log10(0.7425);

lvlEst = [20*log10(vtg100/1000)+Sens; 20*log10(vtg050/1000)+Sens; 20*log10(vtg010/1000)+Sens];

figure;
subplot(3,1,1)
plot( lvlEst(1,:),'bo--','LineWidth',1), hold on, grid on
plot( lvlEst(2,:),'rx--','LineWidth',1)
ylim([101.5 106.5])
xlim([0.8 7.2])
ylabel('SPL [dB]')
ha = gca;
set(ha(end),'XTickLabel',f);
title('(a)')

subplot(3,1,2)
plot( lvlEst(3,:),'bo--','LineWidth',1), hold on, grid on
plot( lvlEst(4,:),'rx--','LineWidth',1)
ylim([101.5 106.5]-6)
xlim([0.8 7.2])
ylabel('SPL [dB]')
ha(end+1) = gca;
set(ha(end),'XTickLabel',f);
title('(b)')

subplot(3,1,3)
plot( lvlEst(5,:),'bo--','LineWidth',1), hold on, grid on
plot( lvlEst(6,:),'rx--','LineWidth',1)
ylim([101.5 106.5]-20)
xlim([0.8 7.2])
ylabel('SPL [dB]')
ha(end+1) = gca;
set(ha(end),'XTickLabel',f);
xlabel('Frequency [Hz]')
title('(c)')

plotOpts.I_FontSize = 10;
h = Figure2paperfigureT(gcf,3,1,plotOpts);

ha(end+1) = gca;
set(ha(end),'XTickLabel',f);

fname = '~/Documenten/Documenten-TUe/01-Text/05-Doc-TUe/lx2015-06-12-update/Figures/freq-response-HD265-new.eps';
Saveas(gcf,fname);
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

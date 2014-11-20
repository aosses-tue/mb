function r20141114_FS_model
% function r20141114_FS_model
%
% 1. Description:
%       Code use to generate and to analyse the first version of the 
%       Fluctuation Strength model (see update on 14/11/2014.
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/11/2014
% Last update on: 18/11/2014 % Update this date manually
% Last use on   : 18/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

pathfigures = [Get_TUe_paths('outputs') 'Figures-20141114-FS' delim];
Mkdir(pathfigures);

close all

%%
N_blocks = 1;
[FSAM FSFM hFS afilesFS] = r20141107_fluctuation(N_blocks);

% h = Figure2paperfigure(hFS(1),6,1);
% Saveas(h,[pathfigures 'stim-FS-AM-new']);
% 
% h = Figure2paperfigure(hFS(2),6,1);
% Saveas(h,[pathfigures 'stim-FS-FM-new']);

figure;
f       = [1 2 4 8 16 32];
f_idx   = 1:length(f);

stdFSAM = std(FSAM);
meanFSAM = mean(FSAM)
plot(f_idx,meanFSAM,'b.-'), hold on, grid on
% plot(f_idx,FSAM,'b.-'), hold on, grid on
errorbar(f_idx,meanFSAM,stdFSAM)

set(gca,'XTickLabel',f);
xlim([min(f_idx)-0.5 max(f_idx)+0.5])
xlabel('Modulation Frequency [Hz]')
ylabel('Fluctuation strength [vacil]')
title('Amplitude Modulated test stimuli')
Saveas(gcf,[pathfigures 'FS-AM-new']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

stdFSFM = std(FSFM);
meanFSFM = mean(FSFM)
plot(f_idx,meanFSFM,'r.-'), hold on, grid on
errorbar(f_idx,meanFSFM,stdFSFM,'r')

set(gca,'XTickLabel',f);
xlim([min(f_idx)-0.5 max(f_idx)+0.5])
xlabel('Modulation Frequency [Hz]')
ylabel('Fluctuation strength [vacil]')
title('Frequency Modulated test stimuli')
Saveas(gcf,[pathfigures 'FS-FM-new']);

%%

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

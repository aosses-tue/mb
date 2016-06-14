function r20160425_experiments_WAE
% function r20160425_experiments_WAE
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 25/04/2016
% Last update on: 25/04/2016 
% Last use on   : 25/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

dir = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-8\';
filesConst = {'piano_multi-A-F-intra-inter-constant-EW.apr', 'piano_multi-A-F-intra-inter-constant-AS.apr'};
filesAdapt = {'piano_A-F-intra-inter-adaptive-20-EW.apr','piano_multi-pilot-inter-adaptive-AS.apr'};
dirfigures = ['D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-04-25-update_WAE\Figures\'];
% S1 = Eric
% S2 = Anne

% % Constant stimuli:
% (1,1): procedure02 % GH05 - JBS36
% (1,2): procedure20 % JBS36 - GH05
% (1,3): procedure24 % JBS36 - JBS51-4544
% (1,4): procedure2_t1_t2 % JBS36 t1 - t2
% (1,5): procedure2_t2_t1 % JBS36 t2 - t1
% (1,6): procedure42 % JBS51-4544 - JBS36

i2plot = [4 5 1 2 3 6]; 
procLabels = {'02','20','24','2 t1-t2','2 t2-t1','42'};
procLabels = procLabels(i2plot);

nSubjects = 2;
nPianos = length(i2plot);

scoresConst = nan(nSubjects,length(i2plot));

for i = 1:nSubjects
    fname = [dir filesConst{i}];
    [o1 o2 procID] = quick_constant_multi(fname);
    [procID idxtmp] = sort(procID);
    o1 = o1(idxtmp);
        
    scoresConst(i,:) = o1(i2plot);
    
end

FontSize = 14;
offsetx = 0.15;
figure
plot((1:nPianos)-offsetx, scoresConst(1,:),'r>','LineWidth',2);hold on, grid on
plot((1:nPianos)+offsetx, scoresConst(2,:),'bo','LineWidth',2);
plot(1:nPianos, mean( scoresConst ),'ks','LineWidth',2);
xlim([0.5 nPianos+0.5]);
ylim([0 103]);
Xlabel('Piano pairs',FontSize)
Ylabel('Scores [%]',FontSize)

ha = gca;
set(ha,'XTick',1:nPianos)
set(ha,'XTickLabel',procLabels)
set(ha,'FontSize',FontSize)
legend({'S1','S2','Average'},'Location','SouthEast')
h = gcf;

%%%

N4ave = 8;
Nreve = 12;
reve = [];
thres = [];

for i = 1:nSubjects
    fname = [dir filesAdapt{i}];
    [xx xx xx out] = quick_staircases(fname);
    
    proc = fieldnames(out.staircases);
    
    for k = 1:length(proc)
        figure;
        % Get_mAFC_reversals(out.staircases.(proc{k}));
        [reve(end+1,:) xx reve2plot staircase2plot] = Get_mAFC_reversals(out.staircases.(proc{k}));
        if i == 1
            plot(staircase2plot,'r'); hold on;
            plot(reve2plot, staircase2plot(reve2plot), 'r>','LineWidth',2)
        else
            plot(staircase2plot,'b'); hold on;
            plot(reve2plot, staircase2plot(reve2plot), 'o','LineWidth',2)
        end
                   
        reve2use = reve(end,end-N4ave+1:end);
        thres(end+1,1) = median(reve2use);
        thres(end  ,2) = str2num(proc{k}(end-1:end)); % last two characters
        thres(end  ,3) = i;
        thres(end  ,4) = min(reve2use);
        thres(end  ,5) = max(reve2use);
        Xlabel('Trial number',FontSize)
        Ylabel('Adaptive parameter, SNR [dB]',FontSize)
        h(end+1) = gcf;
        grid on
        Title(sprintf('S%.0f, %s: Threshold = %.2f [dB]',i,proc{k},thres(end,1)),FontSize)
    end
        
end

dataS1 = thres(1,:);
dataS2 = thres(2:end,:);

[id idxtmp] = sort(dataS2(:,2));
dataS2 = dataS2(idxtmp,:); % only first column

ProcID = {'02','20','24','42'};
nProcs = 4;

errorL1 = dataS1(:,1)-dataS1(:,4);
errorU1 = dataS1(:,5)-dataS1(:,1);

errorL2 = dataS2(:,1)-dataS2(:,4);
errorU2 = dataS2(:,5)-dataS2(:,1);

figure;
% plot(2-offsetx,dataS1(1,1),'r>','LineWidth', 2); hold on, grid on
errorbar(2-offsetx,dataS1(1,1),errorL1,errorU1,'r>','LineWidth', 2); hold on, grid on
% plot([1:length(dataS2(:,1))],dataS2(:,1),'bo','LineWidth',2)
errorbar([1:length(dataS2(:,1))],dataS2(:,1),errorL2,errorU2,'bo','LineWidth',2)
legend({'S1','S2','Average'},'Location','NorthEast')
xlim([0.5 nProcs+0.5]);
ylim([-7 6.5])
ha = gca;
set(ha,'XTick',[1:nProcs]);
set(ha,'XTickLabel',ProcID);
set(ha,'FontSize',FontSize);
Xlabel('Piano pairs',FontSize)
Ylabel('Threshold SNR [dB]',FontSize)
h(end+1) = gcf;

for i = 1:length(h)
    Saveas(h(i),[dirfigures 'fig-' num2str(i)],'epsc')
end
disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

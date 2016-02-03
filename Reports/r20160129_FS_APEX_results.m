function r20160129_FS_APEX_results
% function r20160129_FS_APEX_results
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 28/01/2016
% Last update on: 28/01/2016 
% Last use on   : 28/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc

dirC = 'C:\Users\aosses\Desktop\BATWOMAN-workshop\Experiment-data\Modeling the sensation of fluctuation strength\To-present\';
dir_files = [dirC '2015-11-02T09.27.15\Experiment-results-raw-data' delim];

bAM_fmod = 1;
bAM_fc   = 1;
bAM_SPL  = 1;
bAM_d    = 1;

bFM_fmod = 1;

if bAM_fmod
    fm_pooled1 = [];
    fm_pooled2 = [];

    experiment_name = 'AM-fmod';
    test_fmod = [0 0.25 0.5 1 2 4 8 16 32 64 128];
    idxRef = 6;
    
    %%%
           %  1       9 11   15  17
    % Sidx = [1 2 3 4 5 6  7 8    9 10 11];
    Sidx = [2 3 4 7 9 10]; % except S1, S9, S11, S15, S17, S23
    
    f = {'AM-fm_1.apr'; ...
         'AM-fm_3.apr'; ...
         'AM-fm_5.apr'; ...
         'AM-fm_7.apr'; ...
         'AM-fm_9.apr'; ...
         'AM-fm_11.apr'; ...
         'AM-fm_13.apr'; ...
         'AM-fm_15.apr'; ...
         'AM-fm_17.apr'; ...
         'AM-fm_21.apr'; ...
         'AM-fm_23.apr'};
    f = f(Sidx);
    
    for i = 1:length(f)
        fileres = [dir_files f{i}];
        results1   = Get_APEX_slider_results(fileres);

        fm_pooled1 = [fm_pooled1; results1.scores1];
        % fm_pooled1med(i,:) = prctile(results1.scores1,50);
        fm_pooled2 = [fm_pooled2; results1.scores2];
        % fm_pooled2med(i,:) = prctile(results1.scores2,50);
    end
    
    for j = 1:size(fm_pooled1,2)
        tmp = fm_pooled1(:,j);
        idx1(j) = sum(tmp==200);

        tmp = fm_pooled2(:,j);
        idx2(j) = sum(tmp==200);
    end

    tot1 = [prctile(fm_pooled1,75); prctile(fm_pooled1,50); prctile(fm_pooled1,25)];
    var2latex([fm_pooled1; tot1; idx1]);

    tot2 = [prctile(fm_pooled2,75); prctile(fm_pooled2,50); prctile(fm_pooled2,25)]
    var2latex([fm_pooled2; tot2; idx2]);

    [m1, eL1, eU1, perc1] = Prepare_errorbar_perc(fm_pooled1,25,75);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2,25,75);
    factor2normalise = m1(idxRef)/m2(idxRef);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2*factor2normalise,25,75);
            
    % MeansCombined = mean([MeansStd1; MeansStd2*factor2normalise]);
    
    fmod_points = 1:size(fm_pooled1,2);
    
    text_XLabel = 'fmod [Hz]';
    text_YLabel = 'Fluct. Strength [%]';
    
    figure;
    errorbar(fmod_points-0.1,m1,eL1,eU1,'x'), hold on
    errorbar(fmod_points+0.1,m2,eL2,eU2,'rs'),
    % plot(fmod_points,MeansCombined,'ko--','LineWidth',2)
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    legend('Std1','Std2') %,'Combined')
    title(sprintf('Experiment: %s\n. (Factor applied to Std2 = %.4f)',experiment_name,factor2normalise));
    
    %%% Improving the plots:  
    % %              % 9     15
    % % idx = [1 2 3 4 5 6 7 8    9 10 11];
    % 
    % idx = [1 2 3 4 6 7 9 10 11]; % except S9, S15
    % 
    % tot1 = [prctile(fm_pooled1med(idx,:),75); prctile(fm_pooled1med(idx,:),50); prctile(fm_pooled1med(idx,:),25)];
    % var2latex([tot1]);
    % 
    % tot2 = [prctile(fm_pooled2med(idx,:),75); prctile(fm_pooled2med(idx,:),50); prctile(fm_pooled2med(idx,:),25)]
    % var2latex([tot2]);
end
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bAM_fc
    fm_pooled1 = [];
    fm_pooled2 = [];
    idx1 = []; idx2 = [];
    
    experiment_name = 'AM-fc';
    test_fc = [125 250 500 1000 2000 4000 8000];
    idxRef = 4; % 4 = 1000 Hz; 2 = 250 Hz
    
    %%%
         %  1       9 11   15  17 19 21 23 
    Sidx = [1 2 3 4 5 6  7 8    9 10 11 12];
    % Sidx = [2 3 4 7 9 10]; % except S1, S9, S11, S15, S17, S23
    
    f = {'AM-fc_1.apr'; ...
         'AM-fc_3.apr'; ...
         'AM-fc_5.apr'; ...
         'AM-fc_7.apr'; ...
         'AM-fc_9.apr'; ...
         'AM-fc_11.apr'; ...
         'AM-fc_13.apr'; ...
         'AM-fc_15.apr'; ...
         'AM-fc_17.apr'; ...
         'AM-fc_19.apr'; ...
         'AM-fc_21.apr'; ...
         'AM-fc_23.apr'};
    f = f(Sidx);
    
    for i = 1:length(f)
        fileres = [dir_files f{i}];
        results1   = Get_APEX_slider_results(fileres);

        fm_pooled1 = [fm_pooled1; results1.scores1];
        % fm_pooled1med(i,:) = prctile(results1.scores1,50);
        fm_pooled2 = [fm_pooled2; results1.scores2];
        % fm_pooled2med(i,:) = prctile(results1.scores2,50);
    end
    
    for j = 1:size(fm_pooled1,2)
        tmp = fm_pooled1(:,j);
        idx1(j) = sum(tmp==200);

        tmp = fm_pooled2(:,j);
        idx2(j) = sum(tmp==200);
    end

    tot1 = [prctile(fm_pooled1,75); prctile(fm_pooled1,50); prctile(fm_pooled1,25)];
    var2latex([fm_pooled1; tot1; idx1]);

    tot2 = [prctile(fm_pooled2,75); prctile(fm_pooled2,50); prctile(fm_pooled2,25)]
    var2latex([fm_pooled2; tot2; idx2]);

    [m1, eL1, eU1, perc1] = Prepare_errorbar_perc(fm_pooled1,25,75);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2,25,75);
    factor2normalise = m1(idxRef)/m2(idxRef);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2*factor2normalise,25,75);
            
    % MeansCombined = mean([MeansStd1; MeansStd2*factor2normalise]);
    
    fmod_points = 1:size(fm_pooled1,2);
    
    text_XLabel = 'fc [Hz]';
    text_YLabel = 'Fluct. Strength [%]';
    
    figure;
    errorbar(fmod_points-0.1,m1,eL1,eU1,'x'), hold on
    errorbar(fmod_points+0.1,m2,eL2,eU2,'rs'),
    % plot(fmod_points,MeansCombined,'ko--','LineWidth',2)
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fc)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fc)
    legend('Std1','Std2') %,'Combined')
    title(sprintf('Experiment: %s\n. (Factor applied to Std2 = %.4f)',experiment_name,factor2normalise));
    
    disp('')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bAM_SPL
    fm_pooled1 = [];
    fm_pooled2 = [];
    idx1 = []; idx2 = [];
    
    experiment_name = 'AM-SPL';
    test_SPL = [50 60 70 80 90];
    idxRef = 3; % 5: to normalise; 1 = 50 dB; 3 = 70 dB
    
    %%%
         %  1       9 11   15  17 21 23 
    Sidx = [1 2 3 4 5 6  7 8    9 10 11];
    % Sidx = [2 3 4 7 9 10]; % except S1, S9, S11, S15, S17, S23
    
    f = {'AM-SPL_1.apr'; ...
         'AM-SPL_3.apr'; ...
         'AM-SPL_5.apr'; ...
         'AM-SPL_7.apr'; ...
         'AM-SPL_9.apr'; ...
         'AM-SPL_11.apr'; ...
         'AM-SPL_13.apr'; ...
         'AM-SPL_15.apr'; ...
         'AM-SPL_17.apr'; ...
         'AM-SPL_21.apr'; ...
         'AM-SPL_23.apr'};
    f = f(Sidx);
    
    for i = 1:length(f)
        fileres = [dir_files f{i}];
        results1   = Get_APEX_slider_results(fileres);

        fm_pooled1 = [fm_pooled1; results1.scores1];
        % fm_pooled1med(i,:) = prctile(results1.scores1,50);
        fm_pooled2 = [fm_pooled2; results1.scores2];
        % fm_pooled2med(i,:) = prctile(results1.scores2,50);
    end
    
    for j = 1:size(fm_pooled1,2)
        tmp = fm_pooled1(:,j);
        idx1(j) = sum(tmp==200);

        tmp = fm_pooled2(:,j);
        idx2(j) = sum(tmp==200);
    end

    tot1 = [prctile(fm_pooled1,75); prctile(fm_pooled1,50); prctile(fm_pooled1,25)];
    var2latex([fm_pooled1; tot1; idx1]);

    tot2 = [prctile(fm_pooled2,75); prctile(fm_pooled2,50); prctile(fm_pooled2,25)]
    var2latex([fm_pooled2; tot2; idx2]);

    [m1, eL1, eU1, perc1] = Prepare_errorbar_perc(fm_pooled1,25,75);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2,25,75);
    factor2normalise = m1(idxRef)/m2(idxRef);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2*factor2normalise,25,75);
            
    % MeansCombined = mean([MeansStd1; MeansStd2*factor2normalise]);
    
    fmod_points = 1:size(fm_pooled1,2);
    
    text_XLabel = 'SPL [dB]';
    text_YLabel = 'Fluct. Strength [%]';
    
    figure;
    errorbar(fmod_points-0.1,m1,eL1,eU1,'x'), hold on
    errorbar(fmod_points+0.1,m2,eL2,eU2,'rs'),
    % plot(fmod_points,MeansCombined,'ko--','LineWidth',2)
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_SPL)+0.5])
    ha = gca;
    set(ha,'XTick',1:length(test_SPL))
    set(ha,'XTickLabel',test_SPL)
    legend('Std1','Std2') %,'Combined')
    title(sprintf('Experiment: %s\n. (Factor applied to Std2 = %.4f)',experiment_name,factor2normalise));
    
    disp('')
end
%%%

if bAM_d
    fm_pooled1 = [];
    fm_pooled2 = [];
    idx1 = []; idx2 = [];
    
    experiment_name = 'AM-mdepth';
    test_md = [1 2 4 10 20 40];
    idxRef = 3; % 3 = 4 dB; 6 = 40 dB
    
    %%%
         %  1       9 11   15  17 19 21 23 
    Sidx = [1 2 3 4 5 6  7 8    9 10 11 12];
    % Sidx = [2 3 4 7 9 10]; % except S1, S9, S11, S15, S17, S23
    
    f = {'AM-md_1.apr'; ...
         'AM-md_3.apr'; ...
         'AM-md_5.apr'; ...
         'AM-md_7.apr'; ...
         'AM-md_9.apr'; ...
         'AM-md_11.apr'; ...
         'AM-md_13.apr'; ...
         'AM-md_15.apr'; ...
         'AM-md_17.apr'; ...
         'AM-md_19.apr'; ...
         'AM-md_21.apr'; ...
         'AM-md_23.apr'};
    f = f(Sidx);
    
    for i = 1:length(f)
        fileres = [dir_files f{i}];
        results1   = Get_APEX_slider_results(fileres);

        fm_pooled1 = [fm_pooled1; results1.scores1];
        % fm_pooled1med(i,:) = prctile(results1.scores1,50);
        fm_pooled2 = [fm_pooled2; results1.scores2];
        % fm_pooled2med(i,:) = prctile(results1.scores2,50);
    end
    
    for j = 1:size(fm_pooled1,2)
        tmp = fm_pooled1(:,j);
        idx1(j) = sum(tmp==200);

        tmp = fm_pooled2(:,j);
        idx2(j) = sum(tmp==200);
    end

    tot1 = [prctile(fm_pooled1,75); prctile(fm_pooled1,50); prctile(fm_pooled1,25)];
    var2latex([fm_pooled1; tot1; idx1]);

    tot2 = [prctile(fm_pooled2,75); prctile(fm_pooled2,50); prctile(fm_pooled2,25)]
    var2latex([fm_pooled2; tot2; idx2]);

    [m1, eL1, eU1, perc1] = Prepare_errorbar_perc(fm_pooled1,25,75);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2,25,75);
    factor2normalise = m1(idxRef)/m2(idxRef);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2*factor2normalise,25,75);
            
    % MeansCombined = mean([MeansStd1; MeansStd2*factor2normalise]);
    
    fmod_points = 1:size(fm_pooled1,2);
    
    text_XLabel = 'modulation depth [dB]';
    text_YLabel = 'Fluct. Strength [%]';
    
    figure;
    errorbar(fmod_points-0.1,m1,eL1,eU1,'x'), hold on
    errorbar(fmod_points+0.1,m2,eL2,eU2,'rs'),
    % plot(fmod_points,MeansCombined,'ko--','LineWidth',2)
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_md)+0.5])
    ha = gca;
    set(ha,'XTick',1:length(test_md))
    set(ha,'XTickLabel',test_md)
    legend('Std1','Std2') %,'Combined')
    title(sprintf('Experiment: %s\n. (Factor applied to Std2 = %.4f)',experiment_name,factor2normalise));
    
    disp('')
end


if bFM_fmod
    fm_pooled1 = [];
    fm_pooled2 = [];
    idx1 = []; idx2 = [];
    
    experiment_name = 'FM-fmod';
    test_fmod = [0 0.25 0.5 1 2 4 8 16 32 64 128];
    idxRef = 6;
    
    f = {'FM-fm_2.apr'; ...
         'FM-fm_4.apr'; ...
         'FM-fm_6.apr'; ...
         'FM-fm_8.apr'; ...
         'FM-fm_10.apr'; ...
         'FM-fm_12.apr'; ...
         'FM-fm_14.apr'; ...
         'FM-fm_16.apr'; ...
         'FM-fm_18.apr'; ...
         'FM-fm_20.apr'; ...
         'FM-fm_22.apr'; ...
         'FM-fm_24.apr'};

    for i = 1:length(f)
        fileres = [dir_files f{i}];
        results1   = Get_APEX_slider_results(fileres);

        fm_pooled1 = [fm_pooled1; results1.scores1];
        % fm_pooled1med(i,:) = prctile(results1.scores1,50);
        fm_pooled2 = [fm_pooled2; results1.scores2];
        % fm_pooled2med(i,:) = prctile(results1.scores2,50);
    end
    % [r p g] = a3getresults([dir_files f]);

    for j = 1:size(fm_pooled1,2)
        tmp = fm_pooled1(:,j);
        idx1(j) = sum(tmp==200);

        tmp = fm_pooled2(:,j);
        idx2(j) = sum(tmp==200);
    end

    tot1 = [prctile(fm_pooled1,75); prctile(fm_pooled1,50); prctile(fm_pooled1,25)];
    var2latex([fm_pooled1; tot1; idx1]);

    tot2 = [prctile(fm_pooled2,75); prctile(fm_pooled2,50); prctile(fm_pooled2,25)]
    var2latex([fm_pooled2; tot2; idx2]);

    [m1, eL1, eU1, perc1] = Prepare_errorbar_perc(fm_pooled1,25,75);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2,25,75);
    factor2normalise = m1(idxRef)/m2(idxRef);
    [m2, eL2, eU2, perc2] = Prepare_errorbar_perc(fm_pooled2*factor2normalise,25,75);
    
    fmod_points = 1:size(fm_pooled1,2);
    text_XLabel = 'fmod [Hz]';
    text_YLabel = 'Fluct. Strength [%]';
    
    figure;
    errorbar(fmod_points-0.1,m1,eL1,eU1,'x'), hold on
    errorbar(fmod_points+0.1,m2,eL2,eU2,'rs'),
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    legend('Std1','Std2') %,'Combined')
    title(sprintf('Experiment: %s\n. (Factor applied to Std2 = %.4f)',experiment_name,factor2normalise));
    
%     %%%
%           %  1       9     15  17
%     % idx = [1 2 3 4 5 6 7 8    9 10 11];
% 
%     idx = [2 3 4 6 7 9 10 11]; % except S1, S9, S15, S17
% 
%     tot1 = [prctile(fm_pooled1med(idx,:),75); prctile(fm_pooled1med(idx,:),50); prctile(fm_pooled1med(idx,:),25)];
%     var2latex([tot1]);
% 
%     tot2 = [prctile(fm_pooled2med(idx,:),75); prctile(fm_pooled2med(idx,:),50); prctile(fm_pooled2med(idx,:),25)]
%     var2latex([tot2]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

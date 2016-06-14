function r20160605_processing_results_further
% function r20160605_processing_results_further
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 25/05/2016
% Last update on: 25/05/2016 
% Last use on   : 25/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
bDoTriadic = 0; % TODO: check consistency of the MATRIX
bDoStaircase = 1;

dirmain = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Main-ICRA-v3\';
if bDoTriadic
    dir = [dirmain 'All-Triadic' delim];
    file = {'S01__SCRIPT_Triad_Session1_proc.txt', 'S01__SCRIPT_Triad_Session2_proc.txt'; ...
            'S02__SCRIPT_Triad_Session1_proc.txt', 'S02__SCRIPT_Triad_Session2_proc.txt'; ...
            'S03__SCRIPT_Triad_Session1_proc.txt', 'S03__SCRIPT_Triad_Session2_proc.txt'; ...
            'S04__SCRIPT_Triad_Session1_proc.txt', 'S04__SCRIPT_Triad_Session2_proc.txt'; ...
            'S05__SCRIPT_Triad_Session1_proc.txt', 'S05__SCRIPT_Triad_Session2_proc.txt'; ...
            'S06__SCRIPT_Triad_Session1_proc.txt', 'S06__SCRIPT_Triad_Session2_proc.txt'; ... % inconsistent x 2
            'S07__SCRIPT_Triad_Session1_proc.txt', 'S07__SCRIPT_Triad_Session2_proc.txt'; ...
            'S08__SCRIPT_Triad_Session1_proc.txt', 'S08__SCRIPT_Triad_Session2_proc.txt'; ...
            '', ''; ... % 'S09__SCRIPT_Triad_Session1_proc.txt', 'S09__SCRIPT_Triad_Session2_proc.txt'; ...
            'S10__SCRIPT_Triad_Session1_proc.txt', 'S10__SCRIPT_Triad_Session2_proc.txt'; ... % inconsistent x 1
            'S11__SCRIPT_Triad_Session1_proc.txt', 'S11__SCRIPT_Triad_Session2_proc.txt'; ...
            'S12__SCRIPT_Triad_Session1_proc.txt', 'S12__SCRIPT_Triad_Session2_proc.txt'; ...
            'S13__SCRIPT_Triad_Session1_proc.txt', 'S13__SCRIPT_Triad_Session2_proc.txt'; ...
            'S14__SCRIPT_Triad_Session1_proc.txt', 'S14__SCRIPT_Triad_Session2_proc.txt'; ...
            'S15__SCRIPT_Triad_Session1_proc.txt', 'S15__SCRIPT_Triad_Session2_proc.txt'};

    Msim_S01 = il_get_similarity([dir file{1,1}],[dir file{1,2}]);
    Msim_all = Msim_S01;
    
    idxSubjects = [2:8 10:15]; %  1 was already processed
    
    for i = idxSubjects
        if i < 10
            SubjID = sprintf('S0%.0f',i);
        else
            SubjID = sprintf('S%.0f',i);
        end
        exp2evaluate = sprintf('Msim_%s = il_get_similarity([dir file{%.0f,1}],[dir file{%.0f,2}]);',SubjID,i,i);
        if i == 14
            disp('');
        end
        eval(exp2evaluate);
        
        exp2evaluate = sprintf('Msim_all = Msim_all+Msim_%s;',SubjID);
        eval(exp2evaluate);
        
    end

    % Msim_all = Msim_S01+Msim_S02+Msim_S04+Msim_S08;

    N = 7;
    M = 7;
    for i = 1:N
        Msim_idx(i,:) = i*10*ones(1,M);
    end
    for i = 1:M;
        Msim_idx(:,i) = Msim_idx(:,i)+i;
    end
    
    [Mrank_sim idx_sim] = sort(Msim_all(:),'descend');
    Pairs_sim = Msim_idx(idx_sim);
    
    [Pairs_sim Mrank_sim]
    
    disp('')
    var2latex(Msim_all)
    
    D = [];
    for i = 1:size(Msim_all,1)
        D = [D Msim_all(i,i+1:end)];
    end
    
    squareform(D)
    D = D/max(D);
    
    [Y,eigvals] = cmdscale(D);

    dim2account = 2;
    Yred = Y(:,1:dim2account);
    
    % % The best fit would not be used at this moment:
    % [D,Z] = procrustes(X,Yred); % distance and best fit for multidimensional scaling (only 2D)

    Ntreatments = 7;
    xoffset = 0.025;

    figure;
    plot( Yred(:,1),Yred(:,2),'rd' ); grid on
    labels = num2str((1:Ntreatments)');
    % text(X(:,1)+.05,X(:,2),labels,'Color','b');
    text(Yred(:,1)+xoffset,Yred(:,2),labels,'Color','r');
    xlabel('Dimension 1');
    ylabel('Dimension 2');
    % legend({'2D space'},'Location','SE');

    pairs = Get_pairwise_combinations(1,Ntreatments);

    for i = 1:size(pairs,1)
        data1 = Yred(pairs(i,1),:);
        data2 = Yred(pairs(i,2),:);
        leg4plot{i} = sprintf('%.0f%.0f',pairs(i,1),pairs(i,2));
        distance(i,1) = il_euclidean_dist( data1,data2 );
    end

    for i = 1:size(pairs,1)
        data1 = Y(pairs(i,1),1:3);
        data2 = Y(pairs(i,2),1:3);
        % leg4plot{i} = sprintf('%.0f%.0f',pairs(i,1),pairs(i,2));
        distance3D(i,1) = il_euclidean_dist( data1,data2 );
    end

    for i = 1:size(pairs,1)
        data1 = Y(pairs(i,1),1:4);
        data2 = Y(pairs(i,2),1:4);
        % leg4plot{i} = sprintf('%.0f%.0f',pairs(i,1),pairs(i,2));
        distance4D(i,1) = il_euclidean_dist( data1,data2 );
    end
    
    [pairs distance distance3D distance4D]
    disp('')
    % How to get Pearson correlation
end

Results = [];
if bDoStaircase
    dir = [dirmain 'All-AFC' delim];
    file = {'piano_multi-S01-A-adaptive-S01.apr',1; ...
            'piano_multi-S01-B-adaptive-S01.apr',1; ...
            'piano_multi-S01-C-adaptive-S01.apr',1; ...
            'piano_multi-S01-D-adaptive-S01.apr',1; ...
            'piano_multi-S02-A-adaptive-S02.apr',2; ...
            'piano_multi-S02-B-adaptive-S02.apr',2; ...
            'piano_multi-S02-C-adaptive-S02.apr',2; ...
            'piano_multi-S02-D-adaptive-S02.apr',2; ...
            'piano_multi-S03-A-adaptive-S03.apr',3; ...
            'piano_multi-S03-B-adaptive-S03.apr',3; ...
            'piano_multi-S03-C-adaptive-S03.apr',3; ...
            'piano_multi-S03-D-adaptive-S03.apr',3; ...
            'piano_multi-S04-A-adaptive-S04.apr',4; ...
            'piano_multi-S04-B-adaptive-S04.apr',4; ...
            'piano_multi-S04-C-adaptive-S04.apr',4; ...
            'piano_multi-S04-D-adaptive-S04.apr',4; ...
            'piano_multi-S05-A-adaptive-S05.apr',5; ...
            'piano_multi-S05-B-adaptive-S05.apr',5; ...
            'piano_multi-S05-C-adaptive-S05.apr',5; ...
            'piano_multi-S05-D-adaptive-S05.apr',5; ...
            'piano_multi-S06-A-adaptive-S06.apr',6; ...
            'piano_multi-S06-B-adaptive-S06.apr',6; ...
            'piano_multi-S06-C-adaptive-S06.apr',6; ...
            'piano_multi-S06-D-adaptive-S06.apr',6; ...
            'piano_multi-S07-A-adaptive-S07.apr',7; ...
            'piano_multi-S07-B-adaptive-S07.apr',7; ...
            'piano_multi-S07-C-adaptive-S07.apr',7; ...
            'piano_multi-S07-D-adaptive-S07.apr',7; ...
            'piano_multi-S08-A-adaptive-S08.apr',8; ...
            'piano_multi-S08-B-adaptive-S08.apr',8; ...
            'piano_multi-S08-C-adaptive-S08.apr',8; ...
            'piano_multi-S08-D-adaptive-S08.apr',8; ...
            'piano_multi-S09-A-adaptive-S09.apr',9; ...
            'piano_multi-S09-B-adaptive-S09.apr',9; ...
            % 'piano_multi-S09-C-adaptive-S09.apr',9; ... % not measured yet
            % 'piano_multi-S09-D-adaptive-S09.apr',9; ... % not measured yet
%             'piano_multi-S10-A-adaptive-S10.apr',10; ...
%             'piano_multi-S10-B-adaptive-S10.apr',10; ...
%             'piano_multi-S10-C-adaptive-S10.apr',10; ...
%             'piano_multi-S10-D-adaptive-S10.apr',10; ...
            'piano_multi-S11-A-adaptive-S11.apr',11; ...
            'piano_multi-S11-B-adaptive-S11.apr',11; ...
            'piano_multi-S11-C-adaptive-S11.apr',11; ...
            'piano_multi-S11-D-adaptive-S11.apr',11; ...
            'piano_multi-S12-A-adaptive-S12.apr',12; ...
            'piano_multi-S12-B-adaptive-S12.apr',12; ...
            'piano_multi-S12-C-adaptive-S12.apr',12; ...
            'piano_multi-S12-D-adaptive-S12.apr',12; ...
            'piano_multi-S13-A-adaptive-S13.apr',13; ...
            'piano_multi-S13-B-adaptive-S13.apr',13; ...
            'piano_multi-S13-C-adaptive-S13.apr',13; ...
            'piano_multi-S13-D-adaptive-S13.apr',13; ...
            'piano_multi-S14-C-adaptive-S14.apr',14; ...
            'piano_multi-S14-D-adaptive-S14.apr',14; ...
%             'piano_multi-S14-C-adaptive-S14-1.apr',14; ...
%             'piano_multi-S14-D-adaptive-S14-1.apr',14; ...
            'piano_multi-S15-A-adaptive-S15.apr',15}; % ...
            % 'piano_multi-S15-B-adaptive-S15.apr',15}; % this with problems, proc16 with 8 reversals, proc52 with 4
            % 'piano_multi-S15-C-adaptive-S15.apr',15; ...
            % 'piano_multi-S15-D-adaptive-S15.apr',15};
    
    for i = 1:size(file,1)
        opts.mode = 'median';
        opts.N4mean = 8;
        try
            [Thres_tmp xx outs] = quick_staircases([dir file{i,1}],opts);
        catch
            disp('')
        end
        Nproc = length(Thres_tmp);
        for k = 1:Nproc
            idproc(k,1) = str2num( outs.procID{k}(end-1) );
            idproc(k,2) = str2num( outs.procID{k}(end) );
            
            if idproc(k,1) < idproc(k,2)
                idsort(1) = idproc(k,1);
                idsort(2) = idproc(k,2);
                offsetx = 0;
            else
                idsort(2) = idproc(k,1);
                idsort(1) = idproc(k,2);
                offsetx = 1;
            end
            Results(end+1,1:5) = [idproc(k,1)*10+idproc(k,2) idsort(1)*10+idsort(2) Thres_tmp(k) file{i,2} offsetx];
        end
    end
    
    idx2sort = 2; % use 1 to separate 12 and 21.
    [xx idx] = sort(Results(:,idx2sort),'ascend');
    Results_sorted = Results(idx,:);
    
    piano_pair = Results_sorted(:,idx2sort);
    thresholds = Results_sorted(:,3);
    subjects   = Results_sorted(:,4);
    
    piano_unique     = unique(piano_pair);
    piano_unique_idx = 1:length(piano_unique);
    piano_unique_idx = piano_unique_idx(:);
    piano_label      = piano_unique;
    piano_label_idx  = piano_unique_idx;
    
    count_outliers = 0;
    piano_pair_idx = [];
    for i = 1:length(piano_unique)
        idxs = find(piano_pair == piano_unique(i));
        piano_pair_idx(idxs) = i;
        
        idxs2 = find(thresholds(idxs)>20);
        if length(idxs2) ~= 0
            disp('Outlier(s) found...')
            count_outliers = count_outliers + length(idxs2);
            idxs(idxs2) = [];
        end
        
        piano_med(i,1)= median(thresholds(idxs));
        percL = 25;
        percU = 75;
        piano_L(i,1)  = prctile(thresholds(idxs),percL);
        piano_U(i,1)  = prctile(thresholds(idxs),percU);
    end

    offsetx = 0.05;
    offsety = 0.05;
    
    figure; 
    plot(piano_pair_idx        , thresholds,'x','LineWidth',2); grid on, hold on
    text(piano_pair_idx+offsetx, thresholds+offsety,num2str(subjects),'Color','b');
    plot(piano_unique_idx, piano_med,'r>','LineWidth',2);
    ha = gca;
    set(ha,'XTick',piano_label_idx);
    set(ha,'XTickLabel',piano_label);
    ylim([-15 20])
    h = gcf;
    
    figure; 
    errorbar(piano_unique_idx, piano_med, piano_med-piano_L, piano_U-piano_med, 'r>','LineWidth',2);
    grid on
    ha = gca;
    set(ha,'XTick',piano_label_idx);
    set(ha,'XTickLabel',piano_label);
    ylim([-15 20])
    h(end+1) = gcf;
    
    [xx idxSort] = sort(piano_med,'descend');
    
    Trank = piano_unique(idxSort);
    Thres_rank = piano_med(idxSort);

    
    figure; 
    errorbar(piano_unique_idx, piano_med(idxSort), ...
                               piano_med(idxSort)-piano_L(idxSort), ...
                               piano_U(idxSort)-piano_med(idxSort), 'r>','LineWidth',2);
    grid on
    ha = gca;
    set(ha,'XTick',piano_label_idx);
    set(ha,'XTickLabel',piano_label(idxSort));
    ylim([-15 20])
    h(end+1) = gcf;

    
    [Trank Thres_rank piano_L(idxSort) piano_U(idxSort)]
    
    Save_all_figures(h);
    
    disp('')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

function Msim = il_get_similarity(file1,file2);

[triad1 simil1 dissi1] = il_read_txt(file1);
[triad2 simil2 dissi2] = il_read_txt(file2);

triads = [triad1; triad2];
simil  = [simil1; simil2];
dissi  = [dissi1; dissi2];

midsim  = il_get_midsim(simil,dissi);
Msim = il_sim_matrix(triads, simil,midsim);

function [triad simil dissi] = il_read_txt(file);

fid=fopen(file,'r');

count = 1;
idx = 1;
idx_tr = 1:3;
idx_sim_most  = 4:6;
idx_sim_least = 7:9;

while ~feof(fid)
   lin=fgetl(fid);
   if count == 1
        [word1,COUNT,ERRMSG,NEXTINDEX] = sscanf(lin,'%s',9);
   else
        [numbers,COUNT,ERRMSG,NEXTINDEX] = sscanf(lin,'%d',9);
        if length(numbers) == 9
            triad_tmp = transpose( numbers(idx_tr) );
            simil_tmp = transpose( numbers(idx_sim_most) );
            dissi_tmp = transpose( numbers(idx_sim_least) );
            [xx isort] = sort(triad_tmp);
            
            triad(idx,:) = triad_tmp(isort);
            simil(idx,:) = simil_tmp(isort);
            dissi(idx,:) = dissi_tmp(isort);
            
            idx = idx+1;
        end
        
        disp('')
   end
   count = count + 1;
   
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function midsim = il_get_midsim(simil,dissi);

N = size(simil,1);

midsim = nan(size(simil));

for i = 1:N
    
    if simil(i,:) == [1 1 0]
        if dissi(i,:) == [1 0 1]
            midsim(i,:) = [0 1 1];
        else
            midsim(i,:) = [1 0 1];
        end
    end
    
    if simil(i,:) == [1 0 1]
        if dissi(i,:) == [1 1 0]
            midsim(i,:) = [0 1 1];
        else
            midsim(i,:) = [1 1 0];
        end
    end
    
    if simil(i,:) == [0 1 1]
        if dissi(i,:) == [1 1 0]
            midsim(i,:) = [1 0 1];
        else
            midsim(i,:) = [1 1 0];
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Msim = il_sim_matrix(triads,simil,midsim);

N = max(max(triads));
weight_sim    = 2;
weight_midsim = 1;

Msim = zeros(N);

for i = 1:size(triads,1)
    idx = find(simil(i,:)==1);
    row_nr  = triads(i,idx(1));
    cell_nr = triads(i,idx(2));
    Msim(row_nr,cell_nr) = Msim(row_nr,cell_nr) + weight_sim;
    
    idx = find(midsim(i,:)==1);
    row_nr  = triads(i,idx(1));
    cell_nr = triads(i,idx(2));
    Msim(row_nr,cell_nr) = Msim(row_nr,cell_nr) + weight_midsim;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = il_euclidean_dist(x,y)

d = sqrt( sum((x-y).^2) );





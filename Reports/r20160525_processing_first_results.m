function r20160525_processing_first_results
% function r20160525_processing_first_results
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

bDoTriadic = 0;
bDoStaircase = 1;

dirmain = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Main-ICRA-v3\';
if bDoTriadic
    dir = [dirmain 'All-Triadic' delim];
    file = {'S01__SCRIPT_Triad_Session1_proc.txt', 'S01__SCRIPT_Triad_Session2_proc.txt'; ...
            'S02__SCRIPT_Triad_Session1_proc.txt', 'S02__SCRIPT_Triad_Session2_proc.txt'; ...
            'S04__SCRIPT_Triad_Session1_proc.txt', 'S04__SCRIPT_Triad_Session2_proc.txt'; ...
            'S08__SCRIPT_Triad_Session1_proc.txt', 'S08__SCRIPT_Triad_Session2_proc.txt'};

    Msim_S01 = il_get_similarity([dir file{1,1}],[dir file{1,2}]);
    Msim_S02 = il_get_similarity([dir file{2,1}],[dir file{2,2}]);
    Msim_S04 = il_get_similarity([dir file{3,1}],[dir file{3,2}]);
    Msim_S08 = il_get_similarity([dir file{4,1}],[dir file{4,2}]);

    Msim_all = Msim_S01+Msim_S02+Msim_S04+Msim_S08;

    var2latex(Msim_all)
end

if bDoStaircase
    dir = [dirmain 'All-AFC' delim];
    file = {'piano_multi-S01-A-adaptive-S01.apr'};
    
    opts.mode = 'median';
    opts.N4mean = 8;
    [Thres_tmp xx outs] = quick_staircases([dir file{1}],opts);
    
    for k = 1:length(Thres_tmp)
        idproc(k,1) = str2num( outs.procID{k}(end-1) );
        idproc(k,2) = str2num( outs.procID{k}(end) );
    end
    
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

disp('')




function r20141024_update_dau_et_al
% function r20141024_update_dau_et_al
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       r20141024_update_dau_et_al;
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/10/2014
% Last update on: 20/10/2014 % Update this date manually
% Last use on   : 20/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

close all

% Common parameters:
options.dB_SPL_noise    = 77;

bExpIIA1 = 1;
bExpIIA2 = 0;
bExpIIA3 = 0;
bExpIIB1 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bExpIIA1 == 1

Threshold = [];
label_experiment = 'Dau1996b-ExpIIA1';
%% ExpII.A.1

% SPL_test = 74;
SPL_test = 66:4:86;
options.nExperiment = 1;
% options.test_onsets = 95e-3;
options.test_onsets = [95:10:195]*1e-3;
outs85  = demo_dau1996b(options); 

idx     = outs85.out_stim1.idx;

test_onsets = options.test_onsets;

N_conditions = length(test_onsets);
for i = 1:N_conditions
    exp1 = sprintf('ir_stim%.0f = outs85.out_stim%.0f.template(:,idx);',i,i);
    eval(exp1);
end

%%

for i = 1:length(SPL_test)
    options.dB_SPL          = SPL_test(i);
    outstest = demo_dau1996b(options); 

    % Decision
    
    for j = 1:N_conditions
        exp1 = sprintf('template = outstest.out_stim%.0f.template_no_norm(:,idx);',j); % just to see whether it works
        exp2 = sprintf('mue(j,i) = optimaldetector(ir_stim%.0f,template);',j);
        eval(exp1);
        eval(exp2);
    end
    
end

criterion_corr = 6.5;

for j = 1:N_conditions
    exp1 = sprintf('mue%.0f = mue(%.0f,:);',j,j);
    exp2 = sprintf('Threshold(%.0f) = interp1(mue%.0f,SPL_test,criterion_corr);',j,j);
    eval(exp1);
    eval(exp2);    
end

% Saving figures
figure;
plot(test_onsets,Threshold,'o--'), grid on
xlabel('Signal onset relative to masker onset [ms]')
ylabel('Masked threshold [dB]')
title(sprintf('Signal integration experiment. Criterion = %.1f',criterion_corr))

%
filename = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename);

Thres_labels = [test_onsets; Threshold];

filename = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename,'Threshold');
disp(['Variable saved as: ' filename '.mat']);

filename = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename,'Thres_labels');
disp(['Variable saved as: ' filename '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bExpIIA2 == 1

Threshold = [];
label_experiment = 'Dau1996b-ExpIIA2';
%% ExpII.A.2
% SPL_test = 74;

SPL_test = 68:4:84;
options.nExperiment     = 2;
options.dB_SPL          = 85;
% options.test_phases     = 0;
options.test_phases     = 0:2/8:2; % 9 points, afterwords multiplied by pi
test_phases     = options.test_phases;

outs85          = demo_dau1996b(options); 

N_conditions    = length( options.test_phases );
idx             = outs85.out_stim1.idx;

for i = 1:N_conditions
    exp1 = sprintf('ir_stim%.0f = outs85.out_stim%.0f.template(:,idx);',i,i);
    eval(exp1);
end

for i = 1:length(SPL_test)
    
    options.dB_SPL          = SPL_test(i);
    outstest = demo_dau1996b(options); 

    %% Decision
    
    template    = outstest.out_stim1.template_no_norm(:,idx); %
    mue(1,i)      = optimaldetector(ir_stim1,template);
    
    for j = 2:N_conditions
        exp1 = sprintf('template = outstest.out_stim%.0f.template_no_norm(:,idx);',j); % just to see whether it works
        exp2 = sprintf('mue(j,i) = optimaldetector(ir_stim%.0f,template);',j);
        eval(exp1);
        eval(exp2);
    end
end

criterion_corr = 6.5;

for j = 1:N_conditions
    exp1 = sprintf('mue%.0f = mue(%.0f,:);',j,j);
    exp2 = sprintf('Threshold(%.0f) = interp1(mue%.0f,SPL_test,criterion_corr);',j,j);
    eval(exp1);
    eval(exp2);    
end

% Saving figures
figure;
plot(test_phases,Threshold,'o--'), grid on
xlabel('Signal phase [rad x \pi]')
ylabel('Masked threshold [dB]')
title(sprintf('Signal integration experiment. Criterion = %.1f',criterion_corr))

filename = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename);

Thres_labels = [test_phases; Threshold];

filename = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename,'Threshold');
disp(['Variable saved as: ' filename '.mat']);

filename = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename,'Thres_labels');
disp(['Variable saved as: ' filename '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bExpIIA3 == 1
   
Threshold3 = [];
label_experiment = 'Dau1996b-ExpIIA3';
%% ExpII.A.3
SPL_test = 60:2:84;
options.nExperiment     = 3;
options.dB_SPL          = 85;
options.stim_durations  = [10 20 40 80 160 320];

stim_durations  = options.stim_durations;
N_conditions    = length(options.stim_durations);
outs85          = demo_dau1996b(options); 

idx = outs85.out_stim1.idx;
ir_stim1    = outs85.out_stim1.template(:,idx);
ir_stim2    = outs85.out_stim2.template(:,idx);
ir_stim3    = outs85.out_stim3.template(:,idx);
ir_stim4    = outs85.out_stim4.template(:,idx);
ir_stim5    = outs85.out_stim5.template(:,idx);
ir_stim6    = outs85.out_stim6.template(:,idx);

for i = 1:length(SPL_test)
    
    options.dB_SPL          = SPL_test(i);
    outstest = demo_dau1996b(options); 

    %% Decision
    
    template    = outstest.out_stim1.template_no_norm(:,idx); % just to see whether it works
    mue(1,i)      = optimaldetector(ir_stim1,template);
    
    for j = 2:N_conditions
        exp1 = sprintf('template = outstest.out_stim%.0f.template_no_norm(:,idx);',j); % just to see whether it works
        exp2 = sprintf('mue(j,i) = optimaldetector(ir_stim%.0f,template);',j);
        eval(exp1);
        eval(exp2);
    end
    
end

criterion_corr = 6.5;

mue1 = mue(1,:);
mue2 = mue(2,:);
mue3 = mue(3,:);

Threshold3(1) = interp1(mue1,SPL_test,criterion_corr);
Threshold3(2) = interp1(mue2,SPL_test,criterion_corr);
Threshold3(3) = interp1(mue3,SPL_test,criterion_corr);

if N_conditions > 3
    mue4 = mue(4,:);
    mue5 = mue(5,:);
    mue6 = mue(6,:);
    
    Threshold3(4) = interp1(mue4,SPL_test,criterion_corr);
    Threshold3(5) = interp1(mue5,SPL_test,criterion_corr);
    Threshold3(6) = interp1(mue6,SPL_test,criterion_corr);
end

%% Saving figures
figure;
plot(stim_durations,Threshold3,'o--'), grid on
xlabel('Signal duration [ms]')
ylabel('Masked threshold [dB]')
title(sprintf('Signal integration experiment. Criterion = %.1f',criterion_corr))

filename = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename);

Thres_labels = [stim_durations; Threshold3];

filename = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename,'Threshold3');
disp(['Variable saved as: ' filename '.mat']);

filename = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename,'Thres_labels');
disp(['Variable saved as: ' filename '.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
if bExpIIB1 == 1
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

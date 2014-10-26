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
% Last update on: 24/10/2014 % Update this date manually
% Last use on   : 24/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;
Diary(mfilename,bDiary);

close all

% Common parameters:
options.dB_SPL_noise = 77;
criterion_corr = 6.5;

bRunningNoise   = 0;
bDeterThres     = 0;
bIntRepr        = 0;

options.method = 'dau1996a';
bExpIA1 = 1; % Temporal position
bExpIA2 = 1; % Relative phase
bExpIA3 = 1; % Signal integration
bExpIB0 = 1; % Forward masking
bExpIC0 = 1; % Backward masking

h = [];

if bRunningNoise
    nStim = 100;
    method = 'dau1996';
    [m_no, s_no, tstats] = r20141024_running_noise(nStim,method);
    method = 'dau1996a';
    [m_ol, s_ol]         = r20141024_running_noise(nStim,method);
    
    figure;
    plot(   tstats*1000,s_no,'o-'), hold on
	plot(   tstats*1000,s_ol,'rx-','LineWidth',2), grid on
    xlabel('time [ms]')
    ylabel('std [MU]')
    title(sprintf('Standard deviation of BP white noise (Int. representation). N = %.0f',nStim))
    xlim([-15 600])
    legend('overshoot limit','no limitation')
    h(end+1) = gcf;
    Saveas(h(end),'noise-int-repr-std');
    
    
    figure;
    plot(   tstats*1000,m_no,'o-'), hold on
    plot(   tstats*1000,m_ol,'rx-','LineWidth',2), grid on
    xlabel('time [ms]')
    ylabel('mean [MU]')
    title(sprintf('Mean of BP white noise (Int. representation). N = %.0f',nStim))
    xlim([-15 600])
    legend('overshoot limit','no limitation')
    h(end+1) = gcf;
    Saveas(h(end),'noise-int-repr-mean');
    
end

if bDeterThres == 1
    nStimuli = 1;
    r20141024_determine_thres(nStimuli);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bIntRepr == 1
    
    idx = 13; % approx fc = 1000;
    pathaudio = Get_TUe_paths('outputs'); % 'D:\Output\';
    
    lvl_dB = [56 76 85];
    
    filename_N      = [ pathaudio 'dau1996b_expIIB0_noisemasker.wav'];
    filename_Sinf   = [ pathaudio 'dau1996b_expIIB0_stim-10ms-' num2str(lvl_dB(1)) '.wav'];
    filename_Sthr   = [ pathaudio 'dau1996b_expIIB0_stim-10ms-' num2str(lvl_dB(2)) '.wav'];
    filename_Sref   = [ pathaudio 'dau1996b_expIIB0_stim-10ms-' num2str(lvl_dB(3)) '.wav'];
    
    [insig_N  fs]   = Wavread(filename_N);
    [insig_S1 fs] = Wavread(filename_Sinf);
    [insig_S2 fs] = Wavread(filename_Sthr);
	[insig_S3 fs] = Wavread(filename_Sref);
    
    durN = 200e-3;
    durS = 10e-3;
    onsetS = [10 115 200]*1e-3;
    hopsize = 10e-3*fs;
    bufsize = 20e-3*fs;
    
    %% 1.1 Noise, overshoot limiting
    [ir_N  fc t] = Get_internal_representations(insig_N      ,fs,'dau1996a');
    [ir_No     ] = Get_internal_representations(insig_N      ,fs,'dau1996'); % overshoot limit
    
    [statsN tbuf]   = Get_stats_from_audio_excerpt(ir_N(:,idx) , fs, bufsize, hopsize);
    [statsNo]       = Get_stats_from_audio_excerpt(ir_No(:,idx), fs, bufsize, hopsize);
    
    strMean = sprintf('\ntotal no ov=%.2f, ov=%.2f',mean(ir_N(:,idx)),mean(ir_No(:,idx)));
    strStd = sprintf('\ntotal no ov=%.0f, ov=%.0f',std(ir_N(:,idx)),std(ir_No(:,idx)));
    
    ha = [];
    figure;
    subplot(3,1,1)
    plot(t,ir_N(:,idx),t,ir_No(:,idx)); grid on
    title(sprintf('Internal representations, freq. band f_c = %.1f [Hz]',fc(idx)))
    ylabel('[MU]')
    ha(end+1) = gca;
    
    subplot(3,1,2)
    plot(tbuf,statsN.mean,'o-',tbuf,statsNo.mean,'x-'); grid on
    title(sprintf('mean, frame=%.0f,hop=%.0f [ms]%s',bufsize/fs*1000,hopsize/fs*1000,strMean))
    % legend('N, no ol','N, ol = 10')
    ylabel('[MU]')
    ha(end+1) = gca;
    
    subplot(3,1,3)
    plot(tbuf,statsN.std,'o-',tbuf,statsNo.std,'x-'); grid on
    title(sprintf('std, frame=%.0f,hop=%.0f [ms]%s',bufsize/fs*1000,hopsize/fs*1000,strStd));
    % legend('N, no ol','N, ol = 10')
    ylabel('[MU]')
    ha(end+1) = gca;
    linkaxes(ha,'x');
    xlim([0 600e-3])
    
    h = gcf;
    h(end+1) = Figure2paperfigure(h(end));
    close(h(end-1));
    legend('lim=0','lim=10')
    
    Saveas(h(end),[Get_TUe_paths('outputs') 'noise-overshoot-limiting']);
    
    %% 1.2 N + S, overshoot limiting
    ir_cond1 = []; ir_cond1o = [];
    ir_cond2 = []; ir_cond2o = [];
    ir_cond3 = []; ir_cond3o = [];
    for j = 1:3
        
        ha = [];
        haRow1 = [];
        haRow2 = [];
        haRow3 = [];

        figure;
        nFigures = 9;
        nCount = 1;
        
        for i = 1:length(onsetS)
            toPad = Gen_silence(onsetS(i),fs);
            N_added = length(toPad);
            strTag = sprintf('%.0f%.0f',j,i);
            exp1 = sprintf('insig%s = [toPad; insig_S%.0f(1:end-N_added)];',strTag,j); % generate the different onsets
            
            exp2 = sprintf('ir_NpS%.0f  = Get_internal_representations(insig_N+insig%s,fs,''dau1996a'');',i,strTag);
            exp3 = sprintf('ir_NpS%.0fo = Get_internal_representations(insig_N+insig%s,fs,''dau1996'' );',i,strTag);
            
            exp4 = sprintf('IR1=ir_NpS%.0f(:,idx);  IR2=ir_NpS%.0fo(:,idx);',i,i);
            
            exp5 = sprintf('[statsS%.0f  tbuf] = Get_stats_from_audio_excerpt(IR1, fs, bufsize, hopsize);',i);
            exp6 = sprintf('[statsS%.0fo]      = Get_stats_from_audio_excerpt(IR2, fs, bufsize, hopsize);',i);
            
            exp7 = 'm1 = mean(IR1);  m2 = mean(IR2);';
            exp8 = 's1 = std(IR1);   s2 = std(IR2);';
            
            exp9 = sprintf('stats1=statsS%.0f;  stats2=statsS%.0fo;',i,i);
            
            exp10 = sprintf('ir_cond%.0f  = [ir_cond%.0f  IR1];',i,i);
            exp11 = sprintf('ir_cond%.0fo = [ir_cond%.0fo IR2];',i,i);
            
            eval(exp1);
            eval(exp2);
            eval(exp3);
            eval(exp4);
            
            eval(exp5);
            eval(exp6);
            eval(exp7);
            eval(exp8);
            
            eval(exp9);
            eval(exp10);
            eval(exp11);

            bDebug = 0;
            if bDebug
                
                countDebug = 1;
                figure;
                expTemp = sprintf('plot(t,insig%s), grid on;',strTag);
                eval(expTemp)
                title(sprintf('Debug-%.0f,level=%.0f',countDebug,lvl_dB(j)))
                pause()
                countDebug = countDebug+1;
                
                tmp_N = ir_N(:,idx);
                tmp_No = ir_No(:,idx);
                
                expTemp = sprintf('plot(t,tmp_N,t,ir_NpS%.0f(:,idx)), grid on;',i);
                eval(expTemp)
                title(sprintf('Debug-%.0f,level=%.0f',countDebug,lvl_dB(j)))
                pause()
                countDebug = countDebug+1;
                
                expTemp = sprintf('plot(t,tmp_N,t,IR1), grid on;');
                eval(expTemp)
                title(sprintf('Debug-%.0f,level=%.0f',countDebug,lvl_dB(j)))
                pause()
                countDebug = countDebug+1;
                
                expTemp = sprintf('plot(t,tmp_No,t,ir_NpS%.0fo(:,idx)), grid on;',i);
                eval(expTemp)
                title(sprintf('Debug-%.0f,level=%.0f',countDebug,lvl_dB(j)))
                pause()
                
                countDebug = countDebug+1;
                expTemp = sprintf('plot(t,tmp_No,t,IR2), grid on;');
                eval(expTemp)
                title(sprintf('Debug-%.0f,level=%.0f',countDebug,lvl_dB(j)))
                pause()
                countDebug = countDebug+1;
                
                disp('')
                close
                             
            end
            strMean = sprintf('\ntotal no ov=%.2f, ov=%.2f',m1,m2);
            strStd  = sprintf('\ntotal no ov=%.0f, ov=%.0f',s1,s2);

            % figure; plot(t,insig1,t); 
            subplot(nFigures,1,nCount)
            plot(t,IR1,t,IR2); grid on, hold on
            plot([onsetS(i) onsetS(i)],[-500 500],'r--')       
            title(sprintf('Int rep, freq. band f_c = %.1f [Hz]',fc(idx)))
            ylabel('[MU]')
            ha(end+1) = gca;
            haRow1(end+1) = gca;
            nCount = nCount + 1;

            subplot(nFigures,1,nCount)
            plot(tbuf,stats1.mean,'o-',tbuf,stats2.mean,'x-'); grid on, hold on
            plot([onsetS(i) onsetS(i)],[-500 500],'r--')
            title(sprintf('mean, frame=%.0f,hop=%.0f [ms]%s',bufsize/fs*1000,hopsize/fs*1000,strMean))
            ylabel('[MU]')
            ha(end+1) = gca;
            nCount = nCount + 1;
            haRow2(end+1) = gca;

            subplot(nFigures,1,nCount)
            plot(tbuf,stats1.std,'o-',tbuf,stats2.std,'x-'); grid on, hold on
            plot([onsetS(i) onsetS(i)],[-500 500],'r--')
            title(sprintf('std, frame=%.0f,hop=%.0f [ms]%s',bufsize/fs*1000,hopsize/fs*1000,strStd));
            ylabel('[MU]')
            ha(end+1) = gca;
            nCount = nCount + 1;
            haRow3(end+1) = gca;
            
        end
        
        linkaxes(ha,'x');
        xlim([0 600e-3]);

        linkaxes(haRow1,'y');
        set(haRow1(1),'YLim',[-200 700])
        linkaxes(haRow2,'y');
        set(haRow2(1),'YLim',[-200 500])
        linkaxes(haRow3,'y');
        set(haRow3(1),'YLim',[-25 220])

        h = gcf;
        h(end+1) = Figure2paperfigure(h(end),3,3);
        close(h(end-1));
        legend('lim=0','lim=10')

        opts.bScale = 0;
        fname = sprintf('signals-overshoot-limiting-%.0f',j);
        Saveas(h(end),fname,opts);
        
        ir_NpS_per_onset    = [ir_NpS1(:,idx)   ir_NpS2(:,idx)  ir_NpS3(:,idx)];
        ir_NpSo_per_onset   = [ir_NpS1o(:,idx)  ir_NpS2o(:,idx) ir_NpS3o(:,idx)];
        
        switch j
            case 1
                ir_NpSinf   = ir_NpS_per_onset;
                ir_NpSoinf  = ir_NpSo_per_onset;
            case 2
                ir_NpSthr   = ir_NpS_per_onset;
                ir_NpSothr  = ir_NpSo_per_onset;
            case 3
                ir_NpSsup   = ir_NpS_per_onset;
                ir_NpSosup  = ir_NpSo_per_onset;
        end
        
        for i = 1:3
            
            figure;
            subplot(2,1,1)
            plot(t,ir_N(:,idx),t,ir_NpS_per_onset(:,i)); grid on % different 'j' are for different test levels
            title(sprintf('Int. rep., onset=%.0f [ms], test tone level = %.0f [dB]; freq. band f_c = %.1f [Hz]', onsetS(i), lvl_dB(j), fc(idx)))
            legend('N','S+N')

            subplot(2,1,2)
            plot(t,ir_No(:,idx),t,ir_NpSo_per_onset(:,i)); grid on % different 'j' are for different test levels
            title(sprintf('Int. rep. with lim=10, onset=%.0f [ms], test tone level = %.0f [dB]; freq. band f_c = %.1f [Hz]', onsetS(i), lvl_dB(j), fc(idx)))
            legend('N','S+N')
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bExpIA1 == 1

Threshold = [];
label_experiment    = 'Dau1996b-ExpIA1'; disp(label_experiment);
label_figure        = 'Temporal position';

%% ExpI.A.1

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

for j = 1:N_conditions
    exp1 = sprintf('mue%.0f = mue(%.0f,:);',j,j);
    exp2 = sprintf('Threshold(%.0f) = interp1(mue%.0f,SPL_test,criterion_corr);',j,j);
    eval(exp1);
    eval(exp2);    
end

% Saving figures
figure;
plot(test_onsets*1000,Threshold,'o--'), grid on
xlabel('Signal onset relative to masker onset [ms]')
ylabel('Masked threshold [dB]')
title(sprintf('%s. Criterion = %.1f',label_figure,criterion_corr))

%
filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename_Sref);

Thres_labels = [test_onsets; Threshold];

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename_Sref,'Threshold');
disp(['Variable saved as: ' filename_Sref '.mat']);

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename_Sref,'Thres_labels');
disp(['Variable saved as: ' filename_Sref '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bExpIA2 == 1

Threshold = []; 
mue = [];
label_experiment    = 'Dau1996b-ExpIA2'; disp(label_experiment);
label_figure        = 'Relative phase';
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

for j = 1:N_conditions
    exp1 = sprintf('mue%.0f = mue(%.0f,:);',j,j);
    exp2 = sprintf('Threshold(%.0f) = interp1(mue%.0f,SPL_test,criterion_corr,''linear'',''extrap'');',j,j);
    eval(exp1);
    eval(exp2);    
end

% Saving figures
figure;
plot(test_phases,Threshold,'o--'), grid on
xlabel('Signal phase [rad x \pi]')
ylabel('Masked threshold [dB]')
title(sprintf('%s. Criterion = %.1f',label_figure,criterion_corr))

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename_Sref);

Thres_labels = [test_phases; Threshold];

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename_Sref,'Threshold');
disp(['Variable saved as: ' filename_Sref '.mat']);

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename_Sref,'Thres_labels');
disp(['Variable saved as: ' filename_Sref '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bExpIA3 == 1
   
Threshold3 = [];
mue = [];
label_experiment    = 'Dau1996b-ExpIA3'; disp(label_experiment);
label_figure        = 'Signal integration experiment';

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
title(sprintf('%s. Criterion = %.1f', label_figure, criterion_corr))

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename_Sref);

Thres_labels = [stim_durations; Threshold3];

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename_Sref,'Threshold3');
disp(['Variable saved as: ' filename_Sref '.mat']);

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename_Sref,'Thres_labels');
disp(['Variable saved as: ' filename_Sref '.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
if bExpIB0 == 1;

Threshold = [];
mue = [];
label_experiment    = 'Dau1996b-ExpIB0'; disp(label_experiment);
label_figure        = 'Forward masking experiment';

% SPL_test = 74;
SPL_test = 26:10:76;
options.nExperiment = 10;

noise_onset = 0e-3; 
noise_dur   = 200e-3;
        
noise_offset        = noise_onset + noise_dur;
test_dur            = 10e-3;
test_offset_ref     = [-10:10:40]*1e-3;
options.test_onsets = noise_offset-test_dur+test_offset_ref;
% opts.fc_idx         = 1000;

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

for j = 1:N_conditions
    exp1 = sprintf('mue%.0f = mue(%.0f,:);',j,j);
    exp2 = sprintf('Threshold(%.0f) = interp1(mue%.0f,SPL_test,criterion_corr,''linear'',''extrap'');',j,j);
    eval(exp1);
    eval(exp2);    
end

% Saving figures
figure;
plot(test_offset_ref*1000,Threshold,'o--'), grid on
xlabel('Signal offset relative to masker offset [ms]')
ylabel('Masked threshold [dB]')
title(sprintf('%s. Criterion = %.1f',label_figure,criterion_corr))

%
filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename_Sref);

Thres_labels = [test_onsets; Threshold];

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename_Sref,'Threshold');
disp(['Variable saved as: ' filename_Sref '.mat']);

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename_Sref,'Thres_labels');
disp(['Variable saved as: ' filename_Sref '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
if bExpIC0 == 1;

Threshold = [];
mue = [];
label_experiment    = 'Dau1996b-ExpIC0'; disp(label_experiment);
label_figure        = 'Backward masking experiment';

options.nExperiment = 20;

noise_onset = 0e-3; 
noise_dur   = 200e-3;

% SPL_test = 74;
SPL_test = [10 16:10:76];
test_onset_ref  = [-100:30:-40,-20:4:20]*1e-3;
% test_onset_ref  = [-30:10:10]*1e-3;

options.test_onsets = noise_onset+test_onset_ref;
% opts.fc_idx         = 1000;

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



for j = 1:N_conditions
    exp1 = sprintf('mue%.0f = mue(%.0f,:);',j,j);
    exp2 = sprintf('Threshold(%.0f) = interp1(mue%.0f,SPL_test,criterion_corr,''linear'',''extrap'');',j,j);
    eval(exp1);
    eval(exp2);    
end

% Saving figures
figure;
plot(test_onset_ref*1000,Threshold,'o--'), grid on
xlabel('Signal onset relative to masker onset [ms]')
ylabel('Masked threshold [dB]')
title(sprintf('%s. Criterion = %.1f',label_figure,criterion_corr))

%
filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
Saveas(gcf,filename_Sref);

Thres_labels = [test_onsets; Threshold];

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres'];
save(filename_Sref,'Threshold');
disp(['Variable saved as: ' filename_Sref '.mat']);

filename_Sref = [Get_TUe_paths('outputs') label_experiment '_Thres_labels'];
save(filename_Sref,'Thres_labels');
disp(['Variable saved as: ' filename_Sref '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

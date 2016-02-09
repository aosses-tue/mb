function r20160208_update_FS_after_opt(bParts)
% function r20160208_update_FS_after_opt(bParts)
%
% 1. Description:
%       Runs validation process of the fluctuation strength model. It assumes
%       a weighting as a function of the critical band (gzi different from
%       unity) and varies the power of the cross correlation factor and the
%       effective modulation depth to look for the minimised error.
% 
% 2. Stand-alone example:
%       r20160205_update_FS_optimisation;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also r20160205_update_FS_optimisation
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 04/02/2016
% Last update on: 08/02/2016 
% Last use on   : 08/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bParts = [1 1];
end

dirmain = [Get_TUe_data_paths('lx_Text') 'lx2016-02-04-update-validation-FS' delim];
dirout  = [dirmain 'Test-battery' delim];
dirfigs = [dirmain 'Figures' delim];

% % Summary of the list of sounds:
% 30:38 = AM fmod; % 39:45 = AM fc; % 46:50 = AM SPL; % 24:29 = AM mdepth
% 15:23 = FM fmod; % 51:56 = FM fc; % 57:61 = FM SPL; % 62:66 = FM fdev

files = {[dirout '01-ref_fluct.wav']                                        ,1  ; ... % -10 dB
         [dirout '02-fluct_test_bbn_AM_m_100_fmod_004Hz_60_dBSPL.wav']      ,1.8; ...  % 02-03-04
         ['']                                                               ,3; ...
         ['']                                                               ,0.9; ...
         [dirout '05-man001-longer.wav']                                    ,NaN; ...  
         [dirout '06-vrouw002.wav']                                         ,NaN; ... 
         [dirout '07-meas-ac-4-dist-ane-HP.wav']                            ,NaN; ... 
         [dirout '08-model-ac-4-dist-ane.wav']                              ,NaN; ... 
         [dirout '09-w02-gcirc_f80N_N1000_SR40-44100.wav']                  ,NaN; ... 
         [dirout 'track_37_t01_FM_dev700.wav']                              ,2; ... % 10
         [dirout 'track_37_t02_AM_BBN.wav']                                 ,1.6; ...
         [dirout 'track_37_t03_AM_SIN_2kHz.wav']                            ,1.5; ...
         [dirout 'track_37_t04_FM_dev32.wav']                               ,0.3; ... 
         [dirout 'track_37_t05_NBN.wav']                                    ,0.25; ...
         [dirout '15-FM_tone-fm_0.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'] ,0     ; ... % FM fmod
         [dirout '16-FM_tone-fm_0.25-fc_1500-df_700-SPL_70-w_25-fs_44100-N_176400.wav'],0.21; ...
         [dirout '17-FM_tone-fm_0.50-fc_1500-df_700-SPL_70-w_25-fs_44100-N_176400.wav'],0.43; ...
         [dirout '18-FM_tone-fm_1.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_132300.wav'],0.85; ...
         [dirout '19-FM_tone-fm_2.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'] ,1.17; ...
         [dirout '20-FM_tone-fm_4.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'] ,2   ; ...
         [dirout '21-FM_tone-fm_8.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'] ,0.7; ...
         [dirout '22-FM_tone-fm_16.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'],0.27; ...
         [dirout '23-FM_tone-fm_32.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'],0.02; ...
         [dirout '24-AM_tone-fm_4.00-fc_1000-md_1-SPL_70-w_25-fs_44100-N_88200.wav' ],0   ; ... % AM mdepth (only estimates)
         [dirout '25-AM_tone-fm_4.00-fc_1000-md_2-SPL_70-w_25-fs_44100-N_88200.wav' ],0.02; ...
         [dirout '26-AM_tone-fm_4.00-fc_1000-md_4-SPL_70-w_25-fs_44100-N_88200.wav' ],0.14; ...
         [dirout '27-AM_tone-fm_4.00-fc_1000-md_10-SPL_70-w_25-fs_44100-N_88200.wav'],0.35; ...
         [dirout '28-AM_tone-fm_4.00-fc_1000-md_20-SPL_70-w_25-fs_44100-N_88200.wav'],1.24; ...
         [dirout '29-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'],1.37; ...
         [dirout '30-AM_tone-fm_0.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,0; ... % AM fmod
         [dirout '31-AM_tone-fm_0.25-fc_1000-md_40-SPL_70-w_25-fs_44100-N_176400.wav'],0.06; ...
         [dirout '32-AM_tone-fm_0.50-fc_1000-md_40-SPL_70-w_25-fs_44100-N_176400.wav'],0.16; ...
         [dirout '33-AM_tone-fm_1.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_132300.wav'],0.39; ...
         [dirout '34-AM_tone-fm_2.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,0.84;...
         [dirout '35-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,1.25;...
         [dirout '36-AM_tone-fm_8.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,1.3;...
         [dirout '37-AM_tone-fm_16.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'],0.36;...
         [dirout '38-AM_tone-fm_32.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'],0.06;...
         [dirout '39-AM_tone-fm_4.00-fc_125-md_40-SPL_70-w_25-fs_44100-N_88200.wav']  ,1.14;... % AM fc
         [dirout '40-AM_tone-fm_4.00-fc_250-md_40-SPL_70-w_25-fs_44100-N_88200.wav']  ,1.28;...
         [dirout '41-AM_tone-fm_4.00-fc_500-md_40-SPL_70-w_25-fs_44100-N_88200.wav']  ,1.27;...
         [dirout '42-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,1.17;...
         [dirout '43-AM_tone-fm_4.00-fc_2000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,1.38;...
         [dirout '44-AM_tone-fm_4.00-fc_4000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,1.31;...
         [dirout '45-AM_tone-fm_4.00-fc_8000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,0.98;...
         [dirout '46-AM_tone-fm_4.00-fc_1000-md_40-SPL_50-w_25-fs_44100-N_88200.wav'] ,0.68; ... % AM SPL 
         [dirout '47-AM_tone-fm_4.00-fc_1000-md_40-SPL_60-w_25-fs_44100-N_88200.wav'] ,1.1; ...
         [dirout '48-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'] ,1.53; ...
         [dirout '49-AM_tone-fm_4.00-fc_1000-md_40-SPL_80-w_25-fs_44100-N_88200.wav'] ,1.96; ...
         [dirout '50-AM_tone-fm_4.00-fc_1000-md_40-SPL_90-w_25-fs_44100-N_88200.wav'] ,2.37; ...
         [dirout '51-FM_tone-fm_4.00-fc_500-df_200-SPL_70-w_25-fs_44100-N_88200.wav'] ,2.21; ... % FM fc
         [dirout '52-FM_tone-fm_4.00-fc_1000-df_200-SPL_70-w_25-fs_44100-N_88200.wav'],2.17; ...
         [dirout '53-FM_tone-fm_4.00-fc_1500-df_200-SPL_70-w_25-fs_44100-N_88200.wav'],1.89; ...
         [dirout '54-FM_tone-fm_4.00-fc_3000-df_200-SPL_70-w_25-fs_44100-N_88200.wav'],1.15; ...
         [dirout '55-FM_tone-fm_4.00-fc_6000-df_200-SPL_70-w_25-fs_44100-N_88200.wav'],0.41; ...
         [dirout '56-FM_tone-fm_4.00-fc_8000-df_200-SPL_70-w_25-fs_44100-N_88200.wav'],0.14; ...
         [dirout '57-FM_tone-fm_4.00-fc_1500-df_700-SPL_40-w_25-fs_44100-N_88200.wav'],0.82; ... % FM SPL
         [dirout '58-FM_tone-fm_4.00-fc_1500-df_700-SPL_50-w_25-fs_44100-N_88200.wav'],1.08; ...
         [dirout '59-FM_tone-fm_4.00-fc_1500-df_700-SPL_60-w_25-fs_44100-N_88200.wav'],1.32; ...
         [dirout '60-FM_tone-fm_4.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'],1.58; ...
         [dirout '61-FM_tone-fm_4.00-fc_1500-df_700-SPL_80-w_25-fs_44100-N_88200.wav'],1.82; ...
         [dirout '62-FM_tone-fm_4.00-fc_1500-df_16-SPL_70-w_25-fs_44100-N_88200.wav'] ,0.05; ... % FM fdev
         [dirout '63-FM_tone-fm_4.00-fc_1500-df_32-SPL_70-w_25-fs_44100-N_88200.wav'] ,0.2 ; ...
         [dirout '64-FM_tone-fm_4.00-fc_1500-df_100-SPL_70-w_25-fs_44100-N_88200.wav'],0.73; ...
         [dirout '65-FM_tone-fm_4.00-fc_1500-df_300-SPL_70-w_25-fs_44100-N_88200.wav'],1.36; ...
         [dirout '66-FM_tone-fm_4.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav'],1.94; ...
         }; 

Gainfactor = [-10; ... % 01
              -10; ...
              -10+20; ...
              -10-20; ...
              -10; ... % to make speech at 65 dB SPL
              -10; ...
              0; ...
              0; ...
              0; ...
              -10; ...
              -10; ...
              -10; ...
              -10; ...
              -10; ... % 14
               3; ... % 15
               ];

% Model set-up:
dur = 2;
dataset = 99; % 'calibrating'
bShort = 1;
fs = 44100;
N  = dur*44100; % 1 sec at 44100 Hz
model_par = Get_fluctuation_strength_params(N,fs,dataset);

% [pm pk pg] = r20160205_update_FS_optimisation('info_val_2016-02-07-at-08h-22m-37s.mat'); % pg = 0
% [pm pk pg] = r20160205_update_FS_optimisation('info_val_2016-02-08-at-08h-24m-58s.mat'); % pg = 0
% [pm pk pg] = r20160205_update_FS_optimisation('info_val_2016-02-08-at-18h-00m-27s.mat'); % pg = 1
pm = 1.7;
pk = 0.15;
pg = 2;
model_par.p_g = pg;
model_par.p_m = pm;
model_par.p_k = pk;
model_par.cal = 0.25; % only starting point

% 1. Calibration
        
%%% idx = 1
idx = 1;
[insig fs]  = Wavread(files{idx});

insig  = From_dB(Gainfactor(idx)) * insig;
dBSPL  = rmsdb(insig)+100;
Nlength     = round(dur*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

FS2cal = FluctuationStrength_Garcia(insig1,fs,N,model_par);
model_par.cal = model_par.cal / FS2cal;

%%%

bPart1 = 1;
bPart2 = 0;
bCalculate = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart1
    
    if bCalculate

        % idx_validation  = [2:4 15:66]; % 1 is the reference, to calibrate the model
        idx_validation = 1:66;
        numsamples = length(idx_validation);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Starting the calculation:
        sampleNr = 1;
        k = 1;

        irun = 1;
        jrun = 1;

        tic
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. Confirming whether calibration tone gives 1
        %%% idx = 1
        idx = 1;
        [insig fs]  = Wavread(files{idx});
        insig       = From_dB(Gainfactor(idx)) * insig;
        dBSPL(idx)  = rmsdb(insig)+100;
        Nlength     = round(dur*fs);

        Nstart = 1*fs; % after one second
        insig1 = insig(Nstart:Nstart+Nlength-1);

        FStmp = FluctuationStrength_Garcia(insig1,fs,N,model_par); % it has to be 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. Starting validation
        %%% idx = 2

        for counti = 1:length(idx_validation)

            idx         = idx_validation(counti);
            switch idx
                case {2,3,4}
                    [insig fs]  = Wavread(files{2});
                    insig       = From_dB(Gainfactor(idx)) * insig;
                    dBSPL(sampleNr) = rmsdb(insig)+100;
                    Nlength     = round(dur*fs);

                    if bShort
                        Nstart = 1; % after one second
                        insig1 = insig(Nstart:Nstart+Nlength-1);
                    else
                        insig1 = insig;
                    end
                otherwise
                    [insig fs]  = Wavread(files{idx});
                    if idx >= 15
                        insig   = From_dB(Gainfactor(15))*insig;
                    else
                        insig   = From_dB(Gainfactor(idx))*insig;
                    end
                    Nlength     = round(dur*fs);
                    dBSPL(sampleNr)  = rmsdb(insig)+100;
                    % fluct_tmp = FluctuationStrength_Garcia(insig,fs,N,model_par); 

                    if bShort
                        Nstart = 1;
                        insig1 = insig(Nstart:Nstart+Nlength-1);
                    else
                        insig1 = insig;
                    end
                    FStmp = FluctuationStrength_Garcia(insig1,fs,N,model_par);
            end

            FStmp = FluctuationStrength_Garcia(insig1,fs,N,model_par); % it has to be 1
            if length(FStmp) > 1
                FStmp = median(FStmp);
            end
            FS(sampleNr,1) = FStmp;
            FS(sampleNr,2) = FS(sampleNr,1)-files{idx,2};
            FS(sampleNr,3) = files{idx,2};
            sampleNr = sampleNr + 1;
        end

        errortot = FS(:,1)-FS(:,3);

        info_val(k,:) = [irun jrun numsamples sum(errortot)/numsamples sum(abs(errortot))/numsamples];
        k = k + 1;
        disp('')
        toc
    else
        % load('matlab-20160208.mat')
        load('matlab-20160209-at-1534');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:length(idx_validation)
        id = idx_validation(i);
        try
            fprintf('FS for file %.0f = %.4f [vacil] (error = %.4f, SPL = %.1f [dB])\n',id,FS(i,1),FS(i,2),dBSPL(i));
        catch
            disp('');
        end
    end

    hFig = [];

        % 30:38 = AM fmod
        % 39:45 = AM fc
        % 46:50 = AM SPL
        % 24:29 = AM mdepth

    idxif = [30 38; ...
             39 45; ...
             46 50; ...
             24 29];
    labels = {  'fmod [Hz]'  , [0 0.25 0.5 1 2 4 8 16 32]; % fmod
                'fc [Hz]'    , [125 250 500 1000 2000 4000 8000]; %fc
                'SPL [dB]'   , [50 60 70 80 90]; % SPL
                'mdepth [dB]',[1 2 4 10 20 40]; ... % mdepth
             };
    YLabels = 'Fluct. Strength [vacil]';

    % AM as a function of mdepth
    for i = 1:4

        FS_lit = [];
        tmp = exp_fastl2007(i-1);
        switch i
            case 1
                FS_lit = [0 tmp.FS_AM];
            case 4
                FS_lit = interp1(tmp.xtest_AM,tmp.FS_AM,labels{i,2});
            otherwise
                FS_lit = tmp.FS_AM;
        end
        FS_lit = transpose(FS_lit)/100 * tmp.FS100(1);

        idx2plot = idxif(i,1):idxif(i,2);
        labelticks = labels{i,2};
        figure;
        plot(1:length(idx2plot),FS(idx2plot,1),'o-','LineWidth',2); grid on; hold on
        plot(1:length(idx2plot),FS_lit,'rs-.','LineWidth',2)
        ylabel(YLabels);
        xlabel(labels{i,1});
        title('AM')
        set(gca,'XTick',1:length(idx2plot))
        set(gca,'XTickLabel',labelticks);
        hFig(end+1) = gcf;
    end

    % 15:23 = FM fmod
    % 51:56 = FM fc
    % 57:61 = FM SPL
    % 62:66 = FM fdev

    idxif = [15 23; ...
             51 56; ...
             57 61; ...
             62 66];
    labels = {  'fmod [Hz]', [0 0.25 0.5 1 2 4 8 16 32]; % fmod
                'fc [Hz]'  , [500 1000 1500 3000 6000 8000]; %fc
                'SPL [dB]' , [40 50 60 70 80]; % SPL
                'fdev [Hz]', [16 32 100 300 700]; ... % fdev
             };

    for i = 1:4

        tmp = exp_fastl2007(i-1);
        switch i
            case 1
                FS_lit = [0 tmp.FS_FM];
            case {2,4}
                FS_lit = interp1(tmp.xtest_FM,tmp.FS_FM,labels{i,2});
            otherwise
                FS_lit = tmp.FS_FM;
        end
        FS_lit = transpose(FS_lit)/100 * tmp.FS100(2);

        idx2plot = idxif(i,1):idxif(i,2);
        labelticks = labels{i,2};
        figure;
        plot(1:length(idx2plot),FS(idx2plot,1),'o-','LineWidth',2); grid on; hold on
        plot(1:length(idx2plot),FS_lit,'rs-.','LineWidth',2)
        ylabel(YLabels);
        xlabel(labels{i,1});
        title('FM')
        set(gca,'XTick',1:length(idx2plot))
        set(gca,'XTickLabel',labelticks);
        hFig(end+1) = gcf;
    end

    Save_all_figures(hFig,dirfigs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bPart2
    idx2analyse = [4 2 3 57 59 61];
    
    if bCalculate
        sampleNr = 1;
        for idx = idx2analyse
            switch idx
                    case {2,3,4}
                        [insig fs]  = Wavread(files{2});
                        insig       = From_dB(Gainfactor(idx)) * insig;
                        dBSPL(sampleNr) = rmsdb(insig)+100;
                        Nlength     = round(dur*fs);

                        if bShort
                            Nstart = 1; % after one second
                            insig1 = insig(Nstart:Nstart+Nlength-1);
                        else
                            insig1 = insig;
                        end
                    otherwise
                        [insig fs]  = Wavread(files{idx});
                        if idx >= 15
                            insig   = From_dB(Gainfactor(15))*insig;
                        else
                            insig   = From_dB(Gainfactor(idx))*insig;
                        end
                        Nlength     = round(dur*fs);
                        dBSPL(sampleNr)  = rmsdb(insig)+100;
                        % fluct_tmp = FluctuationStrength_Garcia(insig,fs,N,model_par); 

                        if bShort
                            Nstart = 1;
                            insig1 = insig(Nstart:Nstart+Nlength-1);
                        else
                            insig1 = insig;
                        end
                        % FStmp = FluctuationStrength_Garcia(insig1,fs,N,model_par);
                end

                [FStmp fitmp outstmp] = FluctuationStrength_Garcia(insig1,fs,N,model_par); % it has to be 1
                if length(FStmp) > 1
                    FStmp = median(FStmp);
                end
                FS(sampleNr,1) = FStmp;
                FS(sampleNr,2) = FS(sampleNr,1)-files{idx,2};
                FS(sampleNr,3) = files{idx,2};

                Mdept(sampleNr,:) = outstmp.mdept;
                Kp(sampleNr,:) = outstmp.kp;
                Gzi(sampleNr,:) = outstmp.gzi;

                sampleNr = sampleNr + 1;

        end
    else
        % load('matlab-20160208-at-1143.mat');
        load('matlab-20160208-at-1201.mat'); % no limit for mdept
        % load('matlab-20160208-at-1252.mat'); % sign for kk
    end
    
    disp('')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

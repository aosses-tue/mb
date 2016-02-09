function [pm_opt pk_opt pg_opt] = r20160205_update_FS_optimisation(file_load_results);
% function [pm_opt pk_opt pg_opt] = r20160205_update_FS_optimisation(file_load_results);
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
%       See also: r20160208_update_FS_after_opt
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 04/02/2016
% Last update on: 08/02/2016 
% Last use on   : 08/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bCalculate = input('Calculate errors (1 = yes; 0 = no, then saved outputs are loaded): ');
    if bCalculate == 0
        file_load_results = 'info_val_2016-02-08-at-08h-24m-58s.mat';
    end
else
    bCalculate = 0;
end

bReduced = input('Reduced validation (1 = yes; 0 = no): ');

dirout = [Get_TUe_data_paths('lx_Text') 'lx2016-02-04-update-validation-FS' delim 'Test-battery' delim];

if bReduced == 0
    idx_validation  = [2:4 15:66]; % 1 is the reference, to calibrate the model
else
    idx_validation  = [2:4 15 17:22 23 26 28:30 34:37 39 42 46 48 50 51 53 55 57 59 61 65:66];
end
    
% % Pilot 1:
% idx_validation = [2:4]; 

idx_AM_tones = [1 24:50];
idx_AM_1 = 30:38; % fmod
% % Pilot 2:
% idx_validation = idx_AM_1;

idx_FM_1 = 15:23; % fmod
% % Pilot 3:
% idx_validation = idx_FM_1;

idx_AM_noise = [2:4 11 14];
idx_FM_tones = [10 13 15:23 51:66];

numsamples = length(idx_validation);

% 30:38 = AM fmod
% 39:45 = AM fc
% 46:50 = AM SPL
% 24:29 = AM mdepth

% 15:23 = FM fmod
% 51:56 = FM fc
% 57:61 = FM SPL
% 62:66 = FM fdev

files = {[dirout '01-ref_fluct.wav']                                        ,1  ; ... % -10 dB
         [dirout '02-fluct_test_bbn_AM_m_100_fmod_004Hz_60_dBSPL.wav']      ,1.8; ...  % 02-03-04
         ['']                                                               ,3; ...
         ['']                                                               ,0.9; ...
         [dirout ''],NaN; ... 
         [dirout ''],NaN; ...
         [dirout ''],NaN; ...
         [dirout ''],NaN; ...
         [dirout ''],NaN; ...
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

% pg = 0; % :.1:2;
% step = 0.1;

% pg = 1;
% step = 0.15;

step = 0.2;
pg = 0+step:step:2;
pm = [1.35 2]; % 0+step:step:2;
pk = [1.35 2]; % 0+step:step:2;

dataset = 99; % 'calibrating'

bShort = 1;
fs = 44100;
dur = 2;
N  = dur*44100; % 1 sec at 44100 Hz
model_par = Get_fluctuation_strength_params(N,fs,dataset);

if bCalculate == 1
    sampleNr = 1;
    k = 1;
    for irun = 1:length(pm)
        for jrun = 1:length(pk)
            for hrun = 1:length(pg);
                fprintf('Processing combination %.0f, %.0f, %.0f\n',irun,jrun,hrun);
                tic
                % 1. Calibration
                model_par.p_g = pg(hrun);
                model_par.p_m = pm(irun);
                model_par.p_k = pk(jrun);
                model_par.cal = 0.25; % only starting point

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

                idxstart    = sampleNr;
                for counti = 1:length(idx_validation)

                    idx         = idx_validation(counti);
                    switch idx
                        case {2,3,4}
                            [insig fs]  = Wavread(files{2});
                            insig       = From_dB(Gainfactor(idx)) * insig;
                            dBSPL(idx)  = rmsdb(insig)+100;
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
                            dBSPL(idx)  = rmsdb(insig)+100;
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
                    idxend = sampleNr;
                    sampleNr = sampleNr+1;

                end

                errortot = FS(idxstart:idxend,2)-FS(idxstart:idxend,3);

                info_val(k,:) = [irun jrun hrun sum(errortot)/numsamples sum(abs(errortot))/numsamples idxstart idxend numsamples];
                k = k + 1;
                disp('')
                toc
            end
        end
    end

    p = Get_date;
    save(['info_val_' p.date2print],'info_val','pg','pk','pm','FS');

end

if bCalculate == 0
    % file_load_results = ['info_val_2016-02-07-at-08h-22m-37s.mat']; % AM-BBN
    % file_load_results = ['info_val_2016-02-07-at-11h-41m-33s.mat']; % AM-fmod, tones
    % file_load_results = ['info_val_2016-02-07-at-16h-08m-12s.mat']; % FM-fmod, tones
    % file_load_results = ['info_val_2016-02-08-at-08h-24m-58s.mat'];
    
    load(file_load_results);
end

idxmin = find(info_val(:,5) == min(info_val(:,5)));
pm_opt = pm(info_val(idxmin,1));
pk_opt = pk(info_val(idxmin,2));
try
    pg_opt = pg(info_val(idxmin,3));
catch
    pg_opt = pg;
end

if nargout == 0
    [pm_opt pk_opt pg_opt info_val(idxmin,5)]
end

% Resuls pilot with idx_validation = 2:4, pg = 0: info_val_2016-02-07-at-08h-22m-37s.mat
%     pm = 2; pk = 0.1; errorabs = 0.961;
% Resuls pilot with idx_validation = 30:38, pg = 0: info_val_2016-02-07-at-11h-41m-33s.mat
%     pm = 0.5; pk = 2; errorabs = 0.0853;
% Resuls pilot with idx_validation = 15:23, pg = 0: info_val_2016-02-07-at-16h-08m-12s.mat
%     pm = 1.7; pk = 0.4; errorabs = 0.0505;
% Resuls pilot with idx_validation = all, pg = 1: info_val_2016-02-08-at-08h-24m-58s.mat
%     pm = 1.3; pk = 0.1; errorabs = 0.8372;

% Results pilot with idx_validation = all (but selection); pg = 1: info_val_2016-02-08-at-18h-00m-27s.mat
 %     pm = 1.6500; pk = 0.1500; errorabs(0.9907)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function y = r20160205_update_FS
% function y = r20160205_update_FS
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 04/02/2016
% Last update on: 04/02/2016 
% Last use on   : 04/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirout = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-02-04-update-validation-FS\Test-battery\';

idx_validation  = [1:4 15:66];
% 30:38 = AM fmod
% 39:45 = AM fc
% 46:50 = AM SPL
% 24:29 = AM mdepth

% 15:23 = FM fmod
% 51:56 = FM fc
% 57:61 = FM SPL
% 62:66 = FM fdev

idx_real_sounds = [5:9];
idx_others      = [10:14];

bDoRealSounds = 0;
bDoOthers     = 0;

files = {[dirout '01-ref_fluct.wav']; ... % -10 dB
         [dirout '02-fluct_test_bbn_AM_m_100_fmod_004Hz_60_dBSPL.wav']; ...  % 02-03-04
         ['']; ...
         ['']; ...
         [dirout '05-man001-longer.wav']; ... 
         [dirout '06-vrouw001.wav']; ...
         [dirout '07-meas-ac-4-dist-ane-HP.wav']; ...
         [dirout '08-model-ac-4-dist-ane.wav']; ...
         [dirout '09-w02-gcirc_f80N_N1000_SR40-44100.wav']; ...
         [dirout 'track_37_t01_FM_dev700.wav']; ...
         [dirout 'track_37_t02_AM_BBN.wav']; ...
         [dirout 'track_37_t03_AM_SIN_2kHz.wav']; ...
         [dirout 'track_37_t04_FM_dev32.wav']; ...
         [dirout 'track_37_t05_NBN.wav']; ...
         [dirout '15-FM_tone-fm_0.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ... % FM fmod
         [dirout '16-FM_tone-fm_0.25-fc_1500-df_700-SPL_70-w_25-fs_44100-N_176400.wav']; ...
         [dirout '17-FM_tone-fm_0.50-fc_1500-df_700-SPL_70-w_25-fs_44100-N_176400.wav']; ...
         [dirout '18-FM_tone-fm_1.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_132300.wav']; ...
         [dirout '19-FM_tone-fm_2.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '20-FM_tone-fm_4.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '21-FM_tone-fm_8.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '22-FM_tone-fm_16.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '23-FM_tone-fm_32.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '24-AM_tone-fm_4.00-fc_1000-md_1-SPL_70-w_25-fs_44100-N_88200.wav']; ... % AM mdepth
         [dirout '25-AM_tone-fm_4.00-fc_1000-md_2-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '26-AM_tone-fm_4.00-fc_1000-md_4-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '27-AM_tone-fm_4.00-fc_1000-md_10-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '28-AM_tone-fm_4.00-fc_1000-md_20-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '29-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '30-AM_tone-fm_0.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']; ... % AM fmod
         [dirout '31-AM_tone-fm_0.25-fc_1000-md_40-SPL_70-w_25-fs_44100-N_176400.wav']; ...
         [dirout '32-AM_tone-fm_0.50-fc_1000-md_40-SPL_70-w_25-fs_44100-N_176400.wav']; ...
         [dirout '33-AM_tone-fm_1.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_132300.wav']; ...
         [dirout '34-AM_tone-fm_2.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '35-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '36-AM_tone-fm_8.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '37-AM_tone-fm_16.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '38-AM_tone-fm_32.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '39-AM_tone-fm_4.00-fc_125-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];... % AM fc
         [dirout '40-AM_tone-fm_4.00-fc_250-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '41-AM_tone-fm_4.00-fc_500-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '42-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '43-AM_tone-fm_4.00-fc_2000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '44-AM_tone-fm_4.00-fc_4000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '45-AM_tone-fm_4.00-fc_8000-md_40-SPL_70-w_25-fs_44100-N_88200.wav'];...
         [dirout '46-AM_tone-fm_4.00-fc_1000-md_40-SPL_50-w_25-fs_44100-N_88200.wav']; ... % AM SPL 
         [dirout '47-AM_tone-fm_4.00-fc_1000-md_40-SPL_60-w_25-fs_44100-N_88200.wav']; ...
         [dirout '48-AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '49-AM_tone-fm_4.00-fc_1000-md_40-SPL_80-w_25-fs_44100-N_88200.wav']; ...
         [dirout '50-AM_tone-fm_4.00-fc_1000-md_40-SPL_90-w_25-fs_44100-N_88200.wav']; ...
         [dirout '51-FM_tone-fm_4.00-fc_500-df_200-SPL_70-w_25-fs_44100-N_88200.wav']; ... % FM fc
         [dirout '52-FM_tone-fm_4.00-fc_1000-df_200-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '53-FM_tone-fm_4.00-fc_1500-df_200-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '54-FM_tone-fm_4.00-fc_3000-df_200-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '55-FM_tone-fm_4.00-fc_6000-df_200-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '56-FM_tone-fm_4.00-fc_8000-df_200-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '57-FM_tone-fm_4.00-fc_1500-df_700-SPL_40-w_25-fs_44100-N_88200.wav']; ... % FM SPL
         [dirout '58-FM_tone-fm_4.00-fc_1500-df_700-SPL_50-w_25-fs_44100-N_88200.wav']; ...
         [dirout '59-FM_tone-fm_4.00-fc_1500-df_700-SPL_60-w_25-fs_44100-N_88200.wav']; ...
         [dirout '60-FM_tone-fm_4.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '61-FM_tone-fm_4.00-fc_1500-df_700-SPL_80-w_25-fs_44100-N_88200.wav']; ...
         [dirout '62-FM_tone-fm_4.00-fc_1500-df_16-SPL_70-w_25-fs_44100-N_88200.wav']; ... % FM fdev
         [dirout '63-FM_tone-fm_4.00-fc_1500-df_32-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '64-FM_tone-fm_4.00-fc_1500-df_100-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '65-FM_tone-fm_4.00-fc_1500-df_300-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         [dirout '66-FM_tone-fm_4.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
         }; 
     
Gainfactor = [-10; ... % 01
              -10; ...
               20; ...
              -20; ...
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

fs = 44100;
dur = 2;
N  = dur*44100; % 1 sec at 44100 Hz
model_par = Get_fluctuation_strength_params(N,fs);

%%% idx = 1
idx = 1;
[insig fs]  = Wavread(files{idx});

insig       = From_dB(Gainfactor(idx)) * insig;
dBSPL(idx)  = rmsdb(insig)+100;
Nlength     = round(dur*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

Nstart = 1*fs;
insig2 = insig(Nstart:Nstart+Nlength-1);

FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par);
% FS(idx,2) = FluctuationStrength_Garcia(insig2,fs,N,model_par);

%%%%%
idx = 2;
[insig fs]  = Wavread(files{idx});
insig       = From_dB(Gainfactor(idx)) * insig;
dBSPL(idx)  = rmsdb(insig)+100;
Nlength     = round(dur*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

Nstart = 1.8*fs;
insig2 = insig(Nstart:Nstart+Nlength-1);

% fluct_tmp = FluctuationStrength_Garcia(insig1,fs,N,model_par); FS(idx,1) = fluct_tmp(1);
FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par);
% FS(idx,2) = FluctuationStrength_Garcia(insig2,fs,N,model_par);

idx = 3; % 80 dB SPL
insig1 = setdbspl(insig1,80);
dBSPL(idx)  = rmsdb(insig1)+100;
FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par); 

idx = 4; % 40 dB SPL
insig1 = setdbspl(insig1,40);
dBSPL(idx)  = rmsdb(insig1)+100;
FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoRealSounds
    for idx = 5:6
        [insig fs]  = Wavread(files{idx});
        insig       = From_dB(Gainfactor(idx)) * insig;
        dBSPL(idx)  = rmsdb(insig)+100;
        Nlength     = round(dur*fs);

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par);
    end

    for idx = 7:9
        [insig fs]  = Wavread(files{idx});
        lvl         = 70;
        insig       = setdbspl(insig,lvl);
        dBSPL(idx)  = rmsdb(insig)+100;
        Nlength     = round(dur*fs);

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);
        FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par);

        Nstart = 3*44100;
        insig2 = insig(Nstart:Nstart+Nlength-1);
        FS(idx,2) = FluctuationStrength_Garcia(insig2,fs,N,model_par);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDoOthers
    for idx = idx_others
        [insig fs]  = Wavread(files{idx});
        insig       = From_dB(Gainfactor(idx))*insig;
        Nlength     = round(dur*fs);
        dBSPL(idx)  = rmsdb(insig)+100;

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);
        FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par);

        Nstart = 3*44100;
        insig2 = insig(Nstart:Nstart+Nlength-1);
        FS(idx,2) = FluctuationStrength_Garcia(insig2,fs,N,model_par);
    end
end

for idx = [57 61] % 15:66
        
    [insig fs]  = Wavread(files{idx});
    insig       = From_dB(Gainfactor(15))*insig;
    Nlength     = round(dur*fs);
    dBSPL(idx)  = rmsdb(insig)+100;
    % fluct_tmp = FluctuationStrength_Garcia(insig,fs,N,model_par); 
    
    Nstart = 1;
    insig1 = insig(Nstart:Nstart+Nlength-1);
    FS(idx,1) = FluctuationStrength_Garcia(insig1,fs,N,model_par);
    
% Nstart = 1*44100;
% insig2 = insig(Nstart:Nstart+Nlength-1);
% FS(idx,2) = FluctuationStrength_Garcia(insig2,fs,N,model_par);
end

for i = 1:length(files)
    try
        fprintf('FS for file %.0f = %.4f, %.4f [vacil] (SPL = %.1f [dB])\n',i,FS(i,1),FS(i,2),dBSPL(i));
    end
end

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
    idx2plot = idxif(i,1):idxif(i,2);
    labelticks = labels{i,2};
    figure;
    plot(1:length(idx2plot),FS(idx2plot,1),'o-','LineWidth',2); grid on
    ylabel(YLabels);
    xlabel(labels{i,1});
    title('AM')
    set(gca,'XTick',1:length(idx2plot))
    set(gca,'XTickLabel',labelticks);
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
            'fc [Hz]'  , [500 1000 1500 3000 6000]; %fc
            'SPL [dB]' , [40 50 60 70 80]; % SPL
            'fdev [Hz]', [16 32 100 200 700]; ... % fdev
         };

for i = 1:4
    idx2plot = idxif(i,1):idxif(i,2);
    labelticks = labels{i,2};
    figure;
    plot(1:length(idx2plot),FS(idx2plot,1),'o-','LineWidth',2); grid on
    ylabel(YLabels);
    xlabel(labels{i,1});
    title('FM')
    set(gca,'XTick',1:length(idx2plot))
    set(gca,'XTickLabel',labelticks);
end

disp('')


%     % % Imported from function Get_specs:
%     % % handles = {@spec_AM_fmod @AM_fc @AM_md @AM_SPL @FM_fm @FM_fc @FM_df @FM_SPL};
%     % handles = {@spec_AM_fmod};
%     % specs   = cellfun(@(x)x(),handles,'UniformOutput',false);
% 
%     tic;
%     model_AM_fmod = cellfun(@(s)il_process_spec(params,s),fnames);
%     toc;
%     
%     figure;
%     plot( 1:length(tests),model_AM_fmod ); grid on
%     set(gca,'XTick',1:length(tests))
%     set(gca,'XTickLabel',tests);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

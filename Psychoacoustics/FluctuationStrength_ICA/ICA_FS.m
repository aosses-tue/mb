function ICA_FS
% function ICA_FS
%
% 1. Description:
%       dataset = 90 for ERB FB, 99 for Terhardt (Terhardt is sensitive to
%       hearing trhreshold for BBN.
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 05/05/2016
% Last update on: 05/05/2016 
% Last use on   : 05/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc
dir_stim = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-05-ICA\Stimuli\';
bPart1_create_stim = 0;
bPart1_do_all = 0;
bCal = 0; 
bTestTerhardt = 0;

if bPart1_do_all == 1;
    bPart2_AMtones = 0;
    bPart2_AM_bbn = 0;
    bPart3_FM = 0;
    bPart4_real = 1;
else    
    bPart2_AMtones = 0;
    bPart2_AM_bbn = 0;
    bPart3_FM = 1;
    bPart4_real = 0;
end

if bPart1_create_stim == 1
    bSave = 1;
else
    bSave = 0;
end

[filenames FS_theo] = ICA_FS_Create(dir_stim,bSave);

% idx 1. Reference
% idx 2-7. AM tones (fmod = 1 - 32 Hz)
% idx 8-13. AM BBN (fmod = 1 - 32 Hz)
% idx 14-19. FM tones (fmod = 1 - 32 Hz)

if bTestTerhardt
    
    idx = 1; % idx = 10: AM-BBN
    [insig fs]  = Wavread(filenames{idx});
    
    f1 = round(bark2hz([10 12 14]));
    dF = 0.5;
    idxf1 = round(f1/dF);
     
    lvl = 70;
    y1 = setdbspl(Create_sin(f1(1),88200/fs,fs),lvl);
    y2 = setdbspl(Create_sin(f1(2),88200/fs,fs),lvl);
    y3 = setdbspl(Create_sin(f1(3),88200/fs,fs),lvl);
    insig = y1+y2+y3;
    
    k = 13*2;
    
    minmaxf = [min(f1)-100 max(f1)+100];
    
    bIdle = 1;
    % PeripheralHearingSystem_t(transpose(insig),fs,bIdle);
    insig = PeripheralHearingSystem_t(transpose(insig),fs,idle); % since 14/05
    [outsig ei_f_dB freq] = TerhardtExcitationPatterns_v3( insig,fs );
    
    ei_f_in = zeros(1,88200);
    ei_f_in(idxf1) = 70;
    
    S1 = [-27 -27 -27];
    S2 = -24 + 230./f1 + 0.2*lvl;
    z  = round(hz2bark(f1));
    zl = round(hz2bark(f1))-5;
    zu = round(hz2bark(f1))+5;
    lvlm5 = 5*S1 + lvl;
    lvlp5 = 5*S2 + lvl;
    figure;
    for i = 1:3
        plot([zl(i) z(i)],[lvlm5(i) lvl],'r'), hold on
        plot([z(i) zu(i)],[lvl lvlp5(i)],'r')
    end
    
    plot(hz2bark(freq),ei_f_dB(k,:),'x','LineWidth',4), grid on, hold on
    plot(hz2bark(freq),ei_f_in,'r')
    xlim(hz2bark(minmaxf))
    ylim([0 lvl+3])
    
    plot([13 13],[0 lvl+3],'k')
    
    outsig = transpose(outsig);
    outsig_w = transpose( sum(outsig,1) ); 
    
    idxz = 1:47;
    zi = 0.5:0.5:23.5;
    fi = bark2hz(zi);
    idxz2plot = [16 17 19]; % 17 Bark should be here
     
    % K = round(88200*0.5);
    % 
    % for i = 1:length(idxz2plot)
    %     [y ydB(:,i) f] = freqfft2(outsig(:,idxz2plot(i)),K,fs); %,'hanning');
    % end
    
    t = (1:size(outsig,1))/fs;
    offsety = 0.02;
    figure;
    plot(t,outsig(:,idxz2plot(:,1))+offsety), hold on
    plot(t,outsig(:,idxz2plot(:,2))        )
    plot(t,outsig(:,idxz2plot(:,3))-offsety)
    
    % figure
    % plot(hz2bark(f),ydB(:,1:length(idxz2plot)));
    % legend('out','in')
    % ylim([0 25])
    % disp('')
    
end


% Model set-up:
dur = 2;
dataset = 99; % 'calibrating Terhardt'
% dataset = 90; % 'calibrating ERB'

fs = 44100;
N  = dur*44100; % 1 sec at 44100 Hz
model_par = Get_fluctuation_strength_params(N,fs,dataset);

pm = 1.7; %1.7;
pk = 1.7; %1.7; % 0.265; % it has little influence
pg = 0; % not accounted

model_par.p_g = pg;
model_par.p_m = pm;
model_par.p_k = pk;

if bCal  == 1
    % % 1.0 Calibration (idx 1)
    
    idx = 1;
    [insig fs]  = Wavread(filenames{idx});

    % insig  = From_dB(Gainfactor(idx)) * insig;
    dBSPL  = rmsdb(insig)+100;
    Nlength     = round(dur*fs);

    Nstart = 1;
    insig1 = insig(Nstart:Nstart+Nlength-1);
    
    model_par.cal = 0.25; % only starting point
    [FS2cal fi outs] = FluctuationStrength_Garcia(insig1,fs,N,model_par);
    model_par.cal = model_par.cal / FS2cal;
    % Update manually this cal value in Get_fluctuation_strength_params
end

%%% Loading again dataset:
dataset = 0; 
model_par = Get_fluctuation_strength_params(N,fs,dataset);
model_par.p_g = pg;
model_par.p_m = pm;
model_par.p_k = pk;

FS = nan(size(FS_theo));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1. Reference again but with loaded cal value (idx = 1):
idx = 1;
[insig fs]  = Wavread(filenames{idx});

% insig  = From_dB(Gainfactor(idx)) * insig;
dBSPL  = rmsdb(insig)+100;
Nlength     = round(dur*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

FS(idx,1)  = FluctuationStrength_Garcia(insig1,fs,N,model_par);

if bPart2_AMtones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 2. Reference but at 70 dB (idx = 4)
idx = 4;
[insig fs]  = Wavread(filenames{idx});

% insig  = From_dB(Gainfactor(idx)) * insig;
dBSPL  = rmsdb(insig)+100;
Nlength     = round(dur*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

FS(idx,1)  = FluctuationStrength_Garcia(insig1,fs,N,model_par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 3. AM tone, fmod = 16 Hz at 70 dB (idx = 6)
idx = 6;
[insig fs]  = Wavread(filenames{idx});

% insig  = From_dB(Gainfactor(idx)) * insig;
dBSPL  = rmsdb(insig)+100;
Nlength     = round(dur*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

FS(idx,1)  = FluctuationStrength_Garcia(insig1,fs,N,model_par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 4. AM tone, fmod = 1 Hz at 70 dB (idx = 2)
idx = 2;
[insig fs]  = Wavread(filenames{idx});

% insig  = From_dB(Gainfactor(idx)) * insig;
dBSPL  = rmsdb(insig)+100;
Nlength     = round(dur*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

FS(idx,1)  = FluctuationStrength_Garcia(insig1,fs,N,model_par);
end

% % 5. AM BBN, fmod = 4 Hz at 60 dB (idx = 10)
idx = [10 7 11 12];
if bPart2_AM_bbn
    for i = idx
        [insig fs]  = Wavread(filenames{i});

        % insig  = From_dB(Gainfactor(idx)) * insig;
        dBSPL  = rmsdb(insig)+100;
        Nlength     = round(dur*fs);

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS(i,1) fi outs]  = FluctuationStrength_Garcia(insig1,fs,N,model_par);
        
        disp('')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 6. FM tones
idx = [16];
if bPart3_FM
    for i = idx
        [insig fs]  = Wavread(filenames{i});

        % insig  = From_dB(Gainfactor(idx)) * insig;
        dBSPL  = rmsdb(insig)+100;
        Nlength     = round(dur*fs);

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS(i,1) fi outs]  = FluctuationStrength_Garcia(insig1,fs,N,model_par);
        
        disp('')
    end
end

idx = 1:length(filenames);

if bPart1_do_all
    
    % idx = 2:7; % AM tones
    % idx = [14 16:17]; % FM tones
    for i = idx
        [insig fs]  = Wavread(filenames{i});

        % insig  = From_dB(Gainfactor(idx)) * insig;
        dBSPL  = rmsdb(insig)+100;
        Nlength     = round(dur*fs);

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS(i,1) fi outs]  = FluctuationStrength_Garcia(insig1,fs,N,model_par);
        
        figure;
        plot(1:47,outs.mdept,1:47,outs.kp); grid on
        legend('m','k')
        title(num2str(i))
         
        disp('')
    end
end

if bPart4_real
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    files = {'Spr_1_1-mono.wav', 1.11, 24; ... % Ni = 1058400 (24 s), lvl = 56.1 dB, 67.2
             'Spr_1_2-mono.wav', 1.21,  9; ... % Ni =  396900 ( 9 s), lvl = 60.0 dB, 69.4
             'Spr_2-mono.wav'  , 0.38,  0; ... % Ni =       1       , lvl = 63.6 dB, 67.8
             'Tier1-mono.wav'  , 1.77,0.5; ... % Ni =   22050 (0.5 s), lvl = 64.5 dB, 73.4 dB (peak)
             'RR2-mono.wav'    , 0.02,  0; ... % Ni = 1, lvl = 60.1 dB
             'M5-mono.wav'     , 0.56,  1; ... % Ni =   44100 (1 s), lvl = 58.2 dB
             'M6_2-mono'       , 0.21, 26};    % Ni =  1146600 (26 s) lvl = 62.1 dB
    
	idx = size(FS,1);
    for i = 1:length(files)
        [insig fs] = Wavread([dir_stim files{i,1}]);
        Ni = round(fs*files{i,2});
        insig = From_dB(-6)*insig(Ni+1:N+Ni,1);
        SPL(i) = rmsdb(insig)+100;
        [FS(idx+i,1) fi outs] = FluctuationStrength_Garcia(insig, fs, N);
        FS_theo(idx+i,1) = files{i,2};
        
        figure;
        plot(1:47,outs.mdept,1:47,outs.kp); grid on
        legend('m','k')
        title(num2str(idx+i))
    end

end

[FS_theo FS]
sum(abs(FS_theo-FS))

if abs( 100*( FS(idx) - FS_theo(idx) )/FS_theo(idx) ) > 10 % 10\% precision
    error('Adjust LP slope, if FS is less than theo decrease slope, increase it otherwise...')
end

disp('')
%%%

% N = fs*2;
% for i = 1:length(fmod)
%     [filename{i}, insig(:,i)] = AM_sine_savewave(f,dur,fs,fmod(i),mdept,SPL,dir_stim); 
%     [fluct(i) fi outs] = FluctuationStrength_Garcia(insig(:,i), fs, N);
% end
% 
% SPL = 60; % Hz
% for i = 1:length(fmod);
%     insigBBN(:,i) = AM_random_noise(20,16000,SPL,dur,fs,fmod(i),mdept);
%     [fluctBBN(i) fi outs] = FluctuationStrength_Garcia(insigBBN(:,i), fs, N);
% end
% 
% f = 1500; % Hz
% deltaf = 700; % Hz
% SPL = 70; % Hz
% for i = 1:length(fmod);
%     [fname, insigFM(:,i)] = FM_sine_savewave(f,dur,fs,fmod(i),deltaf,SPL,dir_stim);
%     filename{i+1} = fname;
%     [fluctFM(i) fi outs] = FluctuationStrength_Garcia(insigFM(:,i), fs, N);
% end
% 
% fres = 1:length(fmod);
% 
% figure;
% plot(fres,fluct,fres,fluctFM); grid on
% ha = gca;
% 
% set(ha,'XTick',fres);
% set(ha,'XTickLabel',fmod);
% xlabel('Modulation frequency [Hz]');
% ylabel('Fluctuation strength [vacil]');
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

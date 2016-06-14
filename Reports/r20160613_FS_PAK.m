function r20160613_FS_PAK
% function r20160613_FS_PAK
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 13/06/2016
% Last update on: 13/06/2016 
% Last use on   : 13/06/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

bSave = 0;
bPart1_AM   = 1;
bPart2_AM_BBN = 1;
bPart3_FM   = 1;
bPart4_real = 1;
N = 2*44100; % 2 s at 44100 Hz
dur = 4;

Bark_ = 0.5*[1:47];

dir_stim = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-05-ICA\Stimuli-4s\';
Mkdir(dir_stim);

% addpath('D:\MATLAB_git\Psychoacoustics\FluctuationStrength_ICA\')
[filenames FS_theo] = ICA_FS_Create(dir_stim,bSave,dur);
FS   = nan(size(FS_theo));
FSmm = nan(size(FS_theo),2);

dir = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-05-ICA\Stimuli-PAK\';

fname = {   'ICA2016_Speech_01_5_2', 'ICA2016_Speech_01_6_1'; ...
            'ICA2016_Speech_02_7_1', 'ICA2016_Speech_02_8_1'; ...
            'ICA2016_Speech_23_4_1', 'ICA2016_Speech_23_5_1';...
            'ICA2016_Animal_34_2_2', 'ICA2016_Animal_34_3_1'; ... 
            'ICA2016_Music_31_9_1' , 'ICA2016_Music_31_10_1'; ...
            'ICA2016_Music_29_10_2', 'ICA2016_Music_29_11_1'};
   
for i = 1:length(fname)
    file = [dir fname{i,1} '.mat'];

    load(file);

    expr = sprintf('raw = %s;',fname{i,1});
    eval(expr)

    expr = sprintf('time_raw = %s_X;',fname{i,1});
    eval(expr)

    figure;
    subplot(1,2,1)
    plot(time_raw,raw); grid on;
    title(sprintf('%s',fname{i,1}))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file = [dir fname{i,2} '.mat'];

    load(file);

    expr = sprintf('raw3D = %s;',fname{i,2});
    eval(expr)

    [xx idx_max] = max(raw(7:end-7));
    [xx idx_min] = min(raw(7:end-7));
    idx_min = idx_min + 7;
    idx_max = idx_max + 7;
    
    expr = sprintf('freq_raw = %s_X;',fname{i,2});
    eval(expr)

    subplot(1,2,2)
    plot(freq_raw,raw3D( :,idx_max(1) )); hold on;
    plot(freq_raw,raw3D( :,idx_min(1) ),'r--'); grid on;
    xlim([0 24])
    xlabel('[Bark]')
    title(sprintf('%s',fname{i,2}))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1. Reference again but with loaded cal value (idx = 1):
idx = 1;
[insig fs]  = Wavread(filenames{idx});
dBSPL  = rmsdb(insig)+100;

Nlength     = round(dur*fs);
% Nlength = min(Nlength,length(insig));
if Nlength > length(insig)
    insig = [insig; insig; insig; insig];
end

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

[FS_ fi]  = FluctuationStrength_Garcia(insig1,fs,N);
FS(idx,1) = median(FS_);
[xx idx_max] = max(FS_);
[xx idx_min] = min(FS_);
FSmm(idx,1:2) = [FS_(idx_min) FS_(idx_max)];

figure;
plot(Bark_,fi(idx_max,:)); grid on, hold on
plot(Bark_,fi(idx_min,:),'r--'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = [2:7];

if bPart1_AM
    for i = 1:length(idx)
        [insig fs]  = Wavread(filenames{idx(i)});
        dBSPL  = rmsdb(insig)+100;

        Nlength     = round(dur*fs);
        % Nlength = min(Nlength,length(insig));
        if Nlength > length(insig)
            insig = [insig; insig; insig; insig];
        end

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS_ fi]  = FluctuationStrength_Garcia(insig1,fs,N);
        FS(idx(i),1) = median(FS_);
        [xx idx_max] = max(FS_);
        [xx idx_min] = min(FS_);
        FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

        figure;
        plot(Bark_,fi(idx_max,:)); grid on, hold on
        plot(Bark_,fi(idx_min,:),'r--'); 
    end
end

idx = [8:13];
if bPart2_AM_BBN
    for i = 1:length(idx)
        [insig fs]  = Wavread(filenames{idx(i)});
        dBSPL  = rmsdb(insig)+100;

        Nlength     = round(dur*fs);
        % Nlength = min(Nlength,length(insig));
        if Nlength > length(insig)
            insig = [insig; insig; insig; insig];
        end

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS_ fi]  = FluctuationStrength_Garcia(insig1,fs,N);
        FS(idx(i),1) = median(FS_);
        [xx idx_max] = max(FS_);
        [xx idx_min] = min(FS_);
        FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

        figure;
        plot(Bark_,fi(idx_max,:)); grid on, hold on
        plot(Bark_,fi(idx_min,:),'r--'); 
    end
end

idx = [14:19];
if bPart3_FM
    for i = 1:length(idx)
        [insig fs]  = Wavread(filenames{idx(i)});
        dBSPL  = rmsdb(insig)+100;

        Nlength     = round(dur*fs);
        % Nlength = min(Nlength,length(insig));
        if Nlength > length(insig)
            insig = [insig; insig; insig; insig];
        end

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS_ fi]  = FluctuationStrength_Garcia(insig1,fs,N);
        FS(idx(i),1) = median(FS_);
        [xx idx_max] = max(FS_);
        [xx idx_min] = min(FS_);
        FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

        figure;
        plot(Bark_,fi(idx_max,:)); grid on, hold on
        plot(Bark_,fi(idx_min,:),'r--'); 
    end
end

idx = idx(end); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart4_real
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    files = {'Spr_1_1-mono.wav', 1.11, 24; ... % Ni = 1058400 (24 s), lvl = 56.1 dB, 67.2
             'Spr_1_2-mono.wav', 1.21,  9; ... % Ni =  396900 ( 9 s), lvl = 60.0 dB, 69.4
             'Spr_2-mono.wav'  , 0.38,  0; ... % Ni =       1       , lvl = 63.6 dB, 67.8
             'Tier1-mono.wav'  , 1.77,0.5; ... % Ni =   22050 (0.5 s), lvl = 64.5 dB, 73.4 dB (peak)
             % 'RR2-mono.wav'    , 0.02,  0; ... % Ni = 1, lvl = 60.1 dB
             'M5-mono.wav'     , 0.56,  1; ... % Ni =   44100 (1 s), lvl = 58.2 dB
             'M6_2-mono'       , 0.21, 26};    % Ni =  1146600 (26 s) lvl = 62.1 dB
    
    for i = 1:length(files)
        [insig fs] = Wavread([dir_stim files{i,1}]);
        Nlength     = round(dur*fs);
        Ni = round(fs*files{i,2});
        insig = From_dB(-6)*insig(Ni+1:Nlength+Ni,1);
        SPL(i) = rmsdb(insig)+100;
        [FS_ fi_ outs] = FluctuationStrength_Garcia(insig, fs, N);
        FS(idx+i,1) = median(FS_);
        [xx idx_max] = max(FS_);
        [xx idx_min] = min(FS_);
        FSmm(idx+i,1:2) = [FS_(idx_min) FS_(idx_max)];

        FS_theo(idx+i,1) = files{i,2};
        
        md    = min(outs.mdept,ones(size(outs.mdept)));
        
        figure;
        plot(Bark_,md,Bark_,outs.kp); grid on, hold on; warning('Arrange this to be max and min')
        % fi_ = 2*outs.cal*(md.^outs.p_m).*(outs.kp.^outs.p_k);
        % FS_ = 0.5*sum(fi_);
        plot(Bark_,2*fi_(idx_max,:),'k-','LineWidth',2); hold on, grid on
        plot(Bark_,2*fi_(idx_min,:),'r--','LineWidth',2); 

        % legend('m','k')
        title(num2str(idx+i))
    end

end

[FS_theo FS FSmm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
 
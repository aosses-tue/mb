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

bSave = 0; % set to 1 to re-generate artificial stimuli
bSaveFigs      = 1;
bDo_AM         = 0;
bDo_realsounds = 1;
bDo_FM         = 0;
bDo_AMBBN      = 1;
bFittingModel = 1;

FontSize = 14;
h = [];
hname = [];
N = 2*44100; % 2 s at 44100 Hz
dur = 4;
durAM = 2; % it does not have any sense to keep it in 4 s, since FS is constant
durFM = 2;

Bark_ = 0.5*[1:47];

dir_main = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-05-ICA\';
dir_stim = [dir_main 'Stimuli\'];
dir_fig  = [dir_main 'Figures-out\'];
Mkdir(dir_stim);
Mkdir(dir_fig);

% addpath('D:\MATLAB_git\Psychoacoustics\FluctuationStrength_ICA\')
[filenames FS_theo] = ICA_FS_Create(dir_stim,bSave,dur);
FS   = nan(size(FS_theo));
FSmm = nan(size(FS_theo),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir = [dir_main 'Stimuli-PAK\'];

files = {'Spr_1_1-mono.wav', 1.11, 24, 'ICA2016_Speech_01_5_2', 'ICA2016_Speech_01_6_1'; ... % Ni = 1058400 (24 s), lvl = 56.1 dB, 67.2
         'Spr_1_2-mono.wav', 1.21,  9, 'ICA2016_Speech_02_7_1', 'ICA2016_Speech_02_8_1'; ... % Ni =  396900 ( 9 s), lvl = 60.0 dB, 69.4
         'Spr_2-mono.wav'  , 0.38,  0, 'ICA2016_Speech_23_4_1', 'ICA2016_Speech_23_5_1'; ... % Ni =       1       , lvl = 63.6 dB, 67.8
         'M6_2-mono'       , 0.21, 26, 'ICA2016_Music_29_10_2', 'ICA2016_Music_29_11_1'; ... % Ni =  1146600 (26 s) lvl = 62.1 dB
         'M5-mono.wav'     , 0.56,  1, 'ICA2016_Music_31_9_1' , 'ICA2016_Music_31_10_1'; ... % Ni =   44100 (1 s), lvl = 58.2 dB
         'Tier1-mono.wav'  , 1.77,0.5, 'ICA2016_Animal_34_2_2', 'ICA2016_Animal_34_3_1'; ... % Ni =   22050 (0.5 s), lvl = 64.5 dB, 73.4 dB (peak)
         'RR2-mono.wav'    , 0.02,  0, '', ''}; ... % Ni = 1, lvl = 60.1 dB

idx = 8:13;
filesAMBBN = { ...
         %'randomnoise-Fc-8010_BW-15980_Fmod-1_Mdept-100_SPL-60.wav',[],[],'',''; ...
         %'randomnoise-Fc-8010_BW-15980_Fmod-2_Mdept-100_SPL-60.wav',[],[],'',''; ...
         'randomnoise-Fc-8010_BW-15980_Fmod-4_Mdept-100_SPL-60.wav',[],[],'ICA2016_AMBBN_fmod04_60dB_18_1','ICA2016_AMBBN_fmod04_60dB_19_1'; ...
         'randomnoise-Fc-8010_BW-15980_Fmod-8_Mdept-100_SPL-60.wav',[],[],'ICA2016_AMBBN_fmod08_60dB_19_2','ICA2016_AMBBN_fmod08_60dB_20_1'}; % ...
         %'randomnoise-Fc-8010_BW-15980_Fmod-16_Mdept-100_SPL-60.wav',[],[],'',''; ...
         %'randomnoise-Fc-8010_BW-15980_Fmod-32_Mdept-100_SPL-60.wav',[],[],'',''};
     
filesFM={ ...
         %'FM-tone-fc-1500_fmod-1_deltaf-700-SPL-70-dB.wav' , [], [], '', ''; ... 
         %'FM-tone-fc-1500_fmod-2_deltaf-700-SPL-70-dB.wav' , [], [], '', ''; ... 
         'FM-tone-fc-1500_fmod-4_deltaf-700-SPL-70-dB.wav' , [], [], 'ICA2016_FMtone_fmod04_70dB_14_2', 'ICA2016_FMtone_fmod04_70dB_15_1'; ... 
         'FM-tone-fc-1500_fmod-8_deltaf-700-SPL-70-dB.wav' , [], [], 'ICA2016_FMtone_fmod08_70dB_16_1', 'ICA2016_FMtone_fmod08_70dB_17_1'}; % ; ... 
         %'FM-tone-fc-1500_fmod-16_deltaf-700-SPL-70-dB.wav', [], [], '' , ''; ... 
         %'FM-tone-fc-1500_fmod-32_deltaf-700-SPL-70-dB.wav', [], [], '', ''};    

allfiles = [files; filesAMBBN; filesFM];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bFittingModel == 0
    if bDo_realsounds
        for i = 1:size(files,1)
            if length(files{i,4})>0
                file = [dir files{i,4} '.mat'];

                load(file);

                expr = sprintf('raw = %s;',files{i,4});
                eval(expr)

                expr = sprintf('time_raw = %s_X;',files{i,4});
                eval(expr)

                figure;
                subplot(1,2,1)
                plot(time_raw,raw); grid on;
                title(name2figname(sprintf('%s',files{i,4})))
                xlabel('Time [s]')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(files{i,5})>0
                file = [dir files{i,5} '.mat'];

                load(file);

                expr = sprintf('raw3D = %s;',files{i,5});
                eval(expr)

                [xx idx_max] = max(raw(7:end-7));
                [xx idx_min] = min(raw(7:end-7));
                idx_min = idx_min + 7-1;
                idx_max = idx_max + 7-1;

                expr = sprintf('freq_raw = %s_X;',files{i,5});
                eval(expr)

                subplot(1,2,2)
                plot(freq_raw,raw3D( :,idx_max(1) )); hold on;
                plot(freq_raw,raw3D( :,idx_min(1) ),'r--'); grid on;
                xlim([0 24])
                xlabel(['Frequency [Bark]'])
                hname{end+1} = name2figname(sprintf('%s',files{i,5}));
                title(hname{end})
                h(end+1) = gcf;
            end

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if bDo_FM
        idx = [14:19];
        for i = 1:size(filesFM,1)
            file = [dir filesFM{i,4} '.mat'];

            load(file);

            expr = sprintf('raw = %s;',filesFM{i,4});
            eval(expr)

            expr = sprintf('time_raw = %s_X;',filesFM{i,4});
            eval(expr)

            figure;
            subplot(1,2,1)
            plot(time_raw,raw); grid on;
            title(name2figname(sprintf('%s',filesFM{i,4})))
            xlabel('Time [s]')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            file = [dir filesFM{i,5} '.mat'];

            load(file);

            expr = sprintf('raw3D = %s;',filesFM{i,5});
            eval(expr)

            [xx idx_max] = max(raw(7:end-7));
            [xx idx_min] = min(raw(7:end-7));
            idx_min = idx_min + 7-1;
            idx_max = idx_max + 7-1;

            expr = sprintf('freq_raw = %s_X;',filesFM{i,5});
            eval(expr)

            subplot(1,2,2)
            plot(freq_raw,raw3D( :,idx_max(1) )); hold on;
            plot(freq_raw,raw3D( :,idx_min(1) ),'r--'); grid on;
            xlim([0 24])
            xlabel(['Frequency [Bark]'])

            hname{end+1} = name2figname(sprintf('%s',filesFM{i,5}));
            title(hname{end})
            h(end+1) = gcf;
        end    

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if bDo_AMBBN
        % idx = [14:19];
        for i = 1:size(filesAMBBN,1)
            file = [dir filesAMBBN{i,4} '.mat'];

            load(file);

            expr = sprintf('raw = %s;',filesAMBBN{i,4});
            eval(expr)

            expr = sprintf('time_raw = %s_X;',filesAMBBN{i,4});
            eval(expr)

            figure;
            subplot(1,2,1)
            plot(time_raw,raw); grid on;
            title(name2figname(sprintf('%s',filesAMBBN{i,4})))
            xlabel('Time [s]')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            file = [dir filesAMBBN{i,5} '.mat'];

            load(file);

            expr = sprintf('raw3D = %s;',filesAMBBN{i,5});
            eval(expr)

            [xx idx_max] = max(raw(7:end-7));
            [xx idx_min] = min(raw(7:end-7));
            idx_min = idx_min + 7-1;
            idx_max = idx_max + 7-1;

            expr = sprintf('freq_raw = %s_X;',filesAMBBN{i,5});
            eval(expr)

            subplot(1,2,2)
            plot(freq_raw,raw3D( :,idx_max(1) )); hold on;
            plot(freq_raw,raw3D( :,idx_min(1) ),'r--'); grid on;
            xlim([0 24])
            xlabel(['Frequency [Bark]'])
            hname{end+1} = name2figname(sprintf('%s',filesAMBBN{i,5}));
            title(hname{end})
            h(end+1) = gcf;

        end    

    end
end

pm_var = 1.7; % [1 1.5 2]; % [0.5 1 1.5 2];
pk_var = 1.7; % [0.5 1 1.5 2];
idx_col = 1;

for id_m = 1:length(pm_var)
    for id_k = 1:length(pk_var)
        if bFittingModel

            h = [];
            suffix = sprintf('-m-%.0f-k-%.0f',pm_var(id_m)*100,pk_var(id_k)*100);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % 1. Reference again but with loaded cal value (idx = 1):
            idx = 1;
            [insig fs]  = Wavread(filenames{idx});
            dBSPL  = rmsdb(insig)+100;

            %%%
            model_par = Get_fluctuation_strength_params(N,fs);
            model_par.p_m = pm_var(id_m);
            model_par.p_k = pk_var(id_k);
            
            shortname = strsplit(filenames{idx},delim);
            shortname = shortname{end};

            Nlength     = round(durAM*fs);

            Nstart = 1;
            insig1 = insig(Nstart:Nstart+Nlength-1);

            [FS_ fi outs]  = FluctuationStrength_Garcia(insig1,fs,N, model_par);
            model_par.cal = model_par.cal / median(FS_);
            % [FS_ fi outs]  = FluctuationStrength_Garcia(insig1,fs,N, model_par);
            
            FS(idx,idx_col) = median(FS_);
            [xx idx_max] = max(FS_);
            [xx idx_min] = min(FS_);
            FSmm(idx,1:2) = [FS_(idx_min) FS_(idx_max)];

            figure;
            subplot(1,2,1)
            plot(outs.t,FS_); grid on
            xlabel('Time [s]')

            subplot(1,2,2)
            plot(Bark_,fi(idx_max,:)); grid on, hold on
            plot(Bark_,fi(idx_min,:),'r--'); 
            xlabel('Frequency [Bark]')
            xlim([0 max(Bark_)])

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            idx = [2:7];
            if bDo_AM

                for i = 1:length(idx)
                    [insig fs]  = Wavread(filenames{idx(i)});
                    dBSPL  = rmsdb(insig)+100;

                    shortname = strsplit(filenames{idx(i)},delim);
                    shortname = [shortname{end} ', p_m=' num2str(model_par.p_m) '; p_k=' num2str(model_par.p_k)];

                    Nlength     = round(durAM*fs);

                    Nstart = 1;
                    insig1 = insig(Nstart:Nstart+Nlength-1);

                    [FS_ fi outs]  = FluctuationStrength_Garcia(insig1,fs,N,model_par);
                    FS(idx(i),idx_col) = median(FS_);
                    [xx idx_max] = max(FS_);
                    [xx idx_min] = min(FS_);
                    FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

                    figure;
                    subplot(1,2,1)
                    plot(outs.t,FS_); grid on, hold on
                    xlabel('Time [s]')
                    title(name2figname(shortname))

                    subplot(1,2,2)
                    plot(Bark_,fi(idx_max,:)); grid on, hold on
                    plot(Bark_,fi(idx_min,:),'r--'); 
                    xlabel('Frequency [Bark]')
                    xlim([0 max(Bark_)])

                    h(end+1) = gcf;
                    hname{end+1} = shortname;
                end

                xoffset = 0.05;
                Xvar = [1:length(idx)];
                XvarLabel = [1 2 4 8 16 32];
                figure;
                plot( Xvar-xoffset,FS(idx,idx_col) ,'s-','LineWidth',2), grid on, hold on;
                plot( Xvar+xoffset,FS_theo(idx) ,'r<-')
                ha = gca;
                set(ha,'XTick',Xvar);
                set(ha,'XTickLabel',XvarLabel);
                set(ha,'FontSize',FontSize);
                Xlabel('f_m_o_d [Hz]',FontSize)
                Ylabel('Fluctuation strength [vacil]',FontSize)
                Title('AM-tones',FontSize)
                legend({'Our model','Literature'},'FontSize',FontSize)
                h(end+1) = gcf;
                hname{end+1} = ['res-AM-tones' suffix];

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            idx = [14:19];
            % idx = 17; % fmod = 8 Hz
            if bDo_FM
                for i = 1:length(idx)
                    [insig fs]  = Wavread(filenames{idx(i)});
                    dBSPL  = rmsdb(insig)+100;

                    shortname = strsplit(filenames{idx(i)},delim);
                    shortname = [shortname{end} ', p_m=' num2str(model_par.p_m) '; p_k=' num2str(model_par.p_k)];

                    Nlength     = round(durFM*fs);
                    % Nlength = min(Nlength,length(insig));
                    if Nlength > length(insig)
                        insig = [insig; insig; insig; insig];
                    end

                    Nstart = 1;
                    insig1 = insig(Nstart:Nstart+Nlength-1);

                    [FS_ fi outs]  = FluctuationStrength_Garcia(insig1,fs,N,model_par);
                    FS(idx(i),idx_col) = median(FS_);
                    [xx idx_max] = max(FS_);
                    [xx idx_min] = min(FS_);
                    FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

                    figure;
                    subplot(1,2,1)
                    plot(outs.t,FS_); grid on, hold on
                    xlabel('Time [s]')
                    title(name2figname(shortname))

                    subplot(1,2,2)
                    plot(Bark_,fi(idx_max,:)); grid on, hold on
                    plot(Bark_,fi(idx_min,:),'r--'); 
                    xlabel('Frequency [Bark]')
                    xlim([0 max(Bark_)])
                    plot(Bark_,outs.mdept,Bark_,outs.kp,Bark_,outs.gzi);

                    h(end+1) = gcf;
                    hname{end+1} = shortname;
                end

                xoffset = 0.05;
                Xvar = [1:length(idx)];
                XvarLabel = [1 2 4 8 16 32];
                figure;
                plot( Xvar-xoffset,FS(idx,idx_col) ,'s-','LineWidth',2), grid on, hold on;
                plot( Xvar+xoffset,FS_theo(idx) ,'r<-')
                ha = gca;
                set(ha,'XTick',Xvar);
                set(ha,'XTickLabel',XvarLabel);
                set(ha,'FontSize',FontSize);
                Xlabel('f_m_o_d [Hz]',FontSize)
                Ylabel('Fluctuation strength [vacil]',FontSize)
                Title('FM-tones',FontSize)
                legend({'Our model','Literature'},'FontSize',FontSize)
                h(end+1) = gcf;
                hname{end+1} = ['res-FM-tones' suffix];

            end

            idx = 11; % [8:13];
            if bDo_AMBBN
                for i = 1:length(idx)
                    [insig fs]= Wavread(filenames{idx(i)});
                    dBSPL     = rmsdb(insig)+100;

                    shortname = strsplit(filenames{idx(i)},delim);
                    shortname = [shortname{end} ', p_m=' num2str(model_par.p_m) '; p_k=' num2str(model_par.p_k)];

                    Nlength     = round(durAM*fs);

                    Nstart = 1;
                    insig1 = insig(Nstart:Nstart+Nlength-1);

                    [FS_ fi outs]  = FluctuationStrength_Garcia(insig1,fs,N,model_par);
                    FS(idx(i),idx_col) = median(FS_);
                    [xx idx_max] = max(FS_);
                    [xx idx_min] = min(FS_);
                    FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

                    if idx(i) == 11
                        figure;
                        plot(Bark_,fi(idx_max,:)/0.5,'k-','LineWidth',2); grid on, hold on
                        % plot(Bark_,fi(idx_min,:),'r--'); 
                        Xlabel('Frequency [Bark]',FontSize)
                        xlim([0 max(Bark_)])
                        Title('AM BBN, f_m_o_d=8 [Hz]',FontSize)
                        fact = 0.05;
                        plot(Bark_,fact*outs.mdept(idx_max,:).^outs.p_m,'m-.','LineWidth',2); 
                        plot(Bark_,fact*outs.kp(idx_max,:).^outs.p_k ,'b--','LineWidth',2);
                        plot(Bark_,fact*outs.gzi(idx_max,:).^outs.p_g,'rx','LineWidth',2);
                        legend({'f_i','m_i','k_i','gz_i'},'Location','NorthWest')
                        Ylabel('Specific fluctuation strength [vacil/Bark]',FontSize)
                        set(gca,'FontSize',FontSize)
                        h(end+1) = gcf;
                        hname{end+1} = 'AM-BBN-case';
                    end

                end
                %%%

                xoffset = 0.05;
                Xvar = [1:length(idx)];
                XvarLabel = [1 2 4 8 16 32];
                figure;
                plot( Xvar-xoffset,FS(idx,idx_col) ,'s-','LineWidth',2), grid on, hold on;
                plot( Xvar+xoffset,FS_theo(idx) ,'r<-')
                ha = gca;
                set(ha,'XTick',Xvar);
                set(ha,'XTickLabel',XvarLabel);
                set(ha,'FontSize',FontSize);
                Xlabel('f_m_o_d [Hz]',FontSize)
                Ylabel('Fluctuation strength [vacil]',FontSize)
                Title('AM BBN',FontSize)
                legend({'Our model','Literature'},'FontSize',FontSize)
                h(end+1) = gcf;
                hname{end+1} = ['res-AMBBN' suffix];
            end

            idxi = 19; % value from last AM-BBN
            idx  = idxi;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if bDo_realsounds

                for i = 6; % 1:length(files)
                    [insig fs]  = Wavread([dir_stim files{i,1}]);
                    Nlength     = round(dur*fs);
                    Ni          = round(fs*files{i,2});
                    insig       = From_dB(-6)*insig(Ni+1:Nlength+Ni,1);
                    SPL(i)      = rmsdb(insig)+100;
                    [FS_ fi_ outs] = FluctuationStrength_Garcia(insig, fs, N,model_par);
                    FS(idx+i,idx_col) = median(FS_);
                    if length(FS_)>1
                        FS_(end) = []; % excludes the last sample
                    end
                    [xx idx_max] = max(FS_);
                    [xx idx_min] = min(FS_);
                    FSmm(idx+i,1:2) = [FS_(idx_min) FS_(idx_max)];

                    shortname = strsplit(files{i,1},delim);
                    shortname = [shortname{end} ', p_m=' num2str(model_par.p_m) '; p_k=' num2str(model_par.p_k)];

                    FS_theo(idx+i,1) = files{i,2};

                    md    = min(outs.mdept,ones(size(outs.mdept)));

                    if i == 6
                        figure;
                        % plot(Bark_,md,Bark_,outs.kp); grid on, hold on; warning('Arrange this to be max and min')
                        % plot(Bark_,2*fi_(idx_max,:),'k-','LineWidth',2); hold on, grid on
                        % plot(Bark_,2*fi_(idx_min,:),'r--','LineWidth',2); 
                        % xlim([0 max(Bark_)])
                        % title(num2str(idx+i))
                        
                        plot(Bark_,fi_(idx_max,:)/0.5,'k-','LineWidth',2); grid on, hold on
                        % plot(Bark_,fi(idx_min,:),'r--'); 
                        Xlabel('Frequency [Bark]',FontSize)
                        xlim([0 max(Bark_)])
                        Title('Duck''s quacking',FontSize)
                        fact = 1;
                        plot(Bark_,fact*outs.mdept(idx_max,:).^outs.p_m,'m-.','LineWidth',2); 
                        plot(Bark_,fact*outs.kp(idx_max,:).^outs.p_k ,'b--','LineWidth',2);
                        plot(Bark_,fact*outs.gzi(idx_max,:).^outs.p_g,'rx','LineWidth',2);
                        legend({'f_i','m_i','k_i','gz_i'},'Location','NorthWest')
                        Ylabel('Specific fluctuation strength [vacil/Bark]',FontSize)
                        set(gca,'FontSize',FontSize)
                        h(end+1) = gcf;
                        hname{end+1} = 'res-Quacks';
                    end
                    
                end

                idx = [idxi+1:idxi+length(files)];

                BarsL = [FS(idx,idx_col)-FSmm(idx,1)];
                BarsU = [FSmm(idx,2)-FS(idx,idx_col)];
                xoffset = 0.11;
                Xvar = [1:length(idx)];
                XvarLabel = [1 2 23 29 31 34 61];
                figure;
                % plot( Xvar-xoffset,FS(idx) ,'s','LineWidth',2), grid on, hold on;
                errorbar( Xvar-xoffset,FS(idx,idx_col),BarsL,BarsU,'s','LineWidth',2), grid on, hold on;
                plot( Xvar+xoffset,FS_theo(idx) ,'r<')
                plot(6-xoffset,FS(idx(6),idx_col),'s','MarkerFaceColor','b')
                plot(6-7*xoffset,FS(idx(6),idx_col),'k*','MarkerSize',12,'LineWidth',2)
                ha = gca;
                set(ha,'XTick',Xvar);
                set(ha,'XTickLabel',XvarLabel);
                set(ha,'FontSize',FontSize);
                Xlabel('Track Nr.',FontSize)
                Ylabel('Fluctuation strength [vacil]',FontSize)
                Title('Everyday sounds + pink noise',FontSize)
                legend({'Our model','Literature'},'FontSize',FontSize)
                h(end+1) = gcf;
                hname{end+1} = ['res-everyday-sounds' suffix];
                ylim([-0.05 2.45]);

            end
        end
        idx_col = idx_col+1;
        
        if bSaveFigs
            for i = 1:length(h)
                Saveas(h(i),[dir_fig hname{i}]);
                opts.format = 'fig';
                Saveas(h(i),[dir_fig hname{i}],opts);
            end
        end
        close all

    end
end

if bSaveFigs
    for i = 1:length(h)
        Saveas(h(i),[dir_fig hname{i}]);
        opts.format = 'fig';
        Saveas(h(i),[dir_fig hname{i}],opts);
    end
end

FS
[FS_theo FS FSmm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
 
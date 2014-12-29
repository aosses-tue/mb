function [FS, afiles] = r20141126_fluctuation(N_blocks)
% function [FS, afiles] = r20141126_fluctuation(N_blocks)
%
% 1. Description:
% 
% 2. Stand-alone example:
%       r20141107_fluctuation;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/11/2014
% Last update on: 25/11/2014 % Update this date manually
% Last use on   : 25/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    N_blocks = 1;
    close all
end

pathaudio_src = [Get_TUe_data_paths('db_Fastl2007')];
pathaudio   = [Get_TUe_paths('outputs') 'Fastl2007_test_20141126' delim];
pathaudio_D = [Get_TUe_paths('outputs') 'Daniel1997_test_20141126' delim];
pathfigures = [Get_TUe_paths('outputs') 'Figures-20141126'      delim];
Mkdir(pathfigures);

h = [];

N           = 8192*4;%*20; % FluctuationStrength_offline_debug, max dim = 8192*15 = 2.78 s
bDebug      = 0;
count_afiles = 1;

bCreate     = 0;

if bCreate
    
    opts.bDoFluct = 1;
    opts.bDoRough = 0;
    opts.dur = (20*N/44100); %  approx. 200e-3;
    opts.bGen_test_tones = 1;
    
    opts.bPsySound = 1;
    opts.bDoRamp = 0; 
    opts.bDoZeroPadding = 0;
    
    outs = Generate_reference_sounds_Zwicker2007(opts);
    % outs2 = Generate_reference_sounds(opts);
    
end

bDoExp0 = 1; % Reference
bDoExp1 = 0; % Fastl2007, Fig.10.1
bDoExp2 = 1; % Fastl2007, Fig.10.2
bDoExp3 = 0; % Fastl2007, Fig.10.3
bDoExp4 = 0; % Fastl2007, Fig.10.4
bDoExp5 = 0; % Fastl2007, Fig.10.5
bDoExp6 = 0; % Fastl2007, Fig.10.6
bDoExp7 = 0; % Fastl2007, Fig.10.7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp0
    
    filenames   = { 'ref_fluct'};
                
    filename = [pathaudio filenames{1} '.wav'];

    [x Fs] = Wavread(filename);

    starti = 1;
    insig = x(starti:starti + N-1);

    out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
    FS0(1) = out{1};

    disp(sprintf('Exp 0: FS=%.3f [vacils]\t test signal: %s\n',out{1},name2figname(filenames{1})));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp1
    % Fastl2007, Fig 10.1 (it considers XX dB SPL)
    
    filenames = {   'track_35_t01', 'track_35_t02', 'track_35_t03'; ...
                    'track_35_t04', 'track_35_t05', 'track_35_t06'; ...
                    'track_35_t07', 'track_35_t08', 'track_35_t09'};
    fmod = [1 4 16];
    titles = {'AM BBN', 'AM Tone', 'FM Tone'};
    
    for k = 1:3
    
        for j = 1:3
            filename = [pathaudio_src filenames{k,j} '.wav'];
            
            [x Fs] = Wavread(filename);

            starti = 1;
            insig = x(starti:starti + N-1);
            
            out = FluctuationStrength_offline_debug(insig, Fs, N, bDebug);
            FS1(k,j) = out{1};
        
            disp(sprintf('Exp 1: FS=%.3f [vacils]\t test signal: %s\n',out{1},name2figname(filenames{k})))
        end
        
        subplot(1,3,k)
        plot(fmod,FS1(k,:),'o--'), grid on
        xlabel('modulation frequency, [Hz]')
        ylabel('Fluctuation Strength [vacil]')
        title( titles{k} )
                
    end
    h(end+1) = gcf;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp2
    
    filenames = {   'track_35_t02',60; ...
                    'track_35_t05',70; ...
                    'track_35_t08',70};
    fmod = 4;
    SPL  = [40:10:80];
    titles = {'AM BBN', 'AM Tone', 'FM Tone'};
    
    figure;
    for k = 1:3
        
        d_dB = SPL - filenames{k,2};
        filename = [pathaudio_src filenames{k,1} '.wav'];
        [x Fs] = Wavread(filename);
        starti = 1;
            
        for j = 1:length(d_dB)
            
            insig = From_dB(d_dB(j))*x(starti:starti + N-1);
           
            out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
            FS2(k,j) = out{1};
        
            disp(sprintf('Exp 2: FS=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames{k})))
        end
        
        subplot(1,3,k)
        plot(SPL,FS2(k,:)), grid on
        if k == 2
            xlabel('Sound Pressure Level [dB]')
        end
        ylabel('FS [vacils]')
        title( titles{k} )
               
    end
    h(end+1) = gcf;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp3
    
    % Fastl2007, Fig. 10.3
    
    mx = [0 6 10 20 40 50 60 70 80 90 100];
    fmod = 4;
    
    for k = 1:length(mx)
        
        filenames = ['fluct_test_bbn_AM_m_' Num2str(mx(k),3) '_fmod_' Num2str(fmod,3) 'Hz_60_dBSPL'];
        filename = [pathaudio filenames '.wav'];
        
        [x Fs] = Wavread(filename);
        starti = 1;

        insig = x(starti:starti + N-1);

        out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
        FS3(k,1) = out{1};

        disp(sprintf('Exp 3.1: FS=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames)))
    
    end
    
    fc = 1000;
    
    figure;
    subplot(1,2,1)
    plot(mx,FS3(:,1)), grid on
    xlabel('Modulation factor [%]')
    ylabel('Fluctuation strength [vacils]')
    
    for k = 1:length(mx)
        filenames   = ['fluct_test_fc_' Num2str(fc,3) '_AM_m_' Num2str(mx(k),3) '_fmod_' Num2str(fmod,3) 'Hz_60_dBSPL'];
        filename    = [pathaudio filenames '.wav'];
        
        [x Fs]      = Wavread(filename);
        starti      = 1;
            
        insig   = x(starti:starti + N-1);

        out     = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
        FS3(k,2) = out{1};

        disp(sprintf('Exp 3.2: FS=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames)))
        
    end
    
    subplot(1,2,2)
    plot(mx,FS3(:,2)), grid on
    h(end+1) = gcf;
    xlabel('Modulation factor [%]')
    ylabel('Fluctuation strength [vacils]')
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp4
    % Fastl2007, Fig. 10.4
    fc      = [125 250 500  1000 2000 4000 8000];
    fmod    = 4;
    lvl     = 70;
    d       = 40; % dB
    
    for k = 1:length(fc)
        
        filenames = ['fluct_test_fc_' Num2str(fc(k),4) '_AM_d_' Num2str(d,2) '_fmod_' Num2str(fmod,3) 'Hz_' num2str(lvl) '_dBSPL'];
        filename = [pathaudio filenames '.wav'];
        
        [x Fs] = Wavread(filename);
        starti = 1;
            
        insig = x(starti:starti + N-1);

        out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
        FS4(k) = out{1};

        disp(sprintf('Exp 4.1: FS=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames)))
        
    end
    figure;
    subplot(1,2,1)
    plot(1:length(fc),FS4(:),'o-'), grid on
    xlabel('Centre frequency')
    title('AM Sin')
    xlim([1-0.5 length(fc)+0.5])
    set(gca,'XTickLabel',fc);
    
    fcfm = [500  1000 2000 4000 8000];
    dev = 200;
    lvl = 70;
    
    for k = 1:length(fcfm)
        
        filenames = ['fluct_test_fc_' Num2str(fcfm(k),4) '_FM_dev_' Num2str(dev,3) '_fmod_' Num2str(fmod,3) 'Hz_' num2str(lvl) '_dBSPL'];
        filename = [pathaudio filenames '.wav'];
        
        [x Fs] = Wavread(filename);
        starti = 1;
            
        insig = x(starti:starti + N-1);

        out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
        FS4_2(k) = out{1};

        disp(sprintf('Exp 4.2: FS=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames)))
        
    end
    
    subplot(1,2,2)
    plot(1:length(fcfm),FS4_2(:),'o-'), grid on
    xlim([1-0.5 length(fcfm)+0.5])
    set(gca,'XTickLabel',fcfm);
    xlabel('Centre frequency')
    title('FM Sin')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp5
    % Fastl2007, Fig.10.5
    fmod    = 4;
    fdev    = [0 16 30 60 120 150 400 600 800 1000];
    fc      = 1500;
    lvl     = 70;
    
    for j = 1:length(fdev)
    
        filenames = ['fluct_test_fc_' Num2str(fc,4) '_FM_dev_' Num2str(fdev(j),3) '_fmod_' Num2str(fmod,3) 'Hz_' num2str(lvl) '_dBSPL'];
        filename = [pathaudio filenames '.wav'];

        [x Fs] = Wavread(filename);

        starti = 1;
        insig = x(starti:starti + N-1);
        out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version

        FS5(j)  = out{1};

        disp(sprintf('Exp 5: FS=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames)))

    end
    
    figure;
    plot(fdev,FS5), grid on
    xlabel('frequency deviation [Hz]')
    ylabel('Fluctuation strength [vacils]')
    h(end+1) = gcf;
    title(sprintf('FM tone, fc=1500 [Hz], fmod=%.0f [Hz], %.0f dB SPL',fmod,lvl))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp6
    % Fastl2007, Fig 10.6 (it considers XX dB SPL)
    
    filenames = {   'track_36_t01', 'track_36_t02', 'track_36_t03'};
    BW = [2 6 50];
    titles = {'NBN'};
    
    figure;
        
    for j = 1:3
        filename = [pathaudio_src filenames{j} '.wav'];

        [x Fs] = Wavread(filename);

        starti = 1;
        insig = x(starti:starti + N-1);

        out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
        FS6(j) = out{1};

        disp(sprintf('Exp 6: FS=%.3f [vacils]\t test signal: %s\n',out{1},name2figname(filenames{j})))
    end

    figure;
    plot(1:3,FS6,'o--'), grid on
    xlabel('Effective modulation frequency [Hz]')
    ylabel('Fluctuation strength [vacils]')
    title( titles{1} )

    xlim([0.5 3.5])
    set(gca,'XTick',1:3);
    set(gca,'XTickLabel',0.64*BW);
    
    h(end+1) = gcf;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp7
    
    filenames   = { 'track_37_t01_FM_dev700',...
                    'track_37_t02_AM_BBN', ...
                    'track_37_t03_AM_SIN_2kHz',...
                    'track_37_t04_FM_dev32',...
                    'track_37_t05_NBN'};
                
    txtLabels   = {'FM SIN 1','AM BBN','AM SIN','FM SIN 2','NBN'};
    
    figure;
        
    for j = 1:5
        filename = [pathaudio_src filenames{j} '.wav'];

        [x Fs] = Wavread(filename);

        starti = 1;
        insig = x(starti:starti + N-1);

        out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
        FS7(j) = out{1};

        disp(sprintf('Exp 6: FS=%.3f [vacils]\t test signal: %s\n',out{1},name2figname(filenames{j})))
    end

    figure;
    plot(1:5,FS7,'o'), grid on
    xlabel('Sound type')
    ylabel('Fluctuation strength [vacils]')

    xlim([0.5 5.5])
    set(gca,'XTick',1:5);
    set(gca,'XTickLabel',txtLabels);
    
    h(end+1) = gcf;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

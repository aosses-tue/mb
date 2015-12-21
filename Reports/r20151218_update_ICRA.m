function r20151218_update_ICRA
% function r20151218_update_ICRA
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 07/12/2015
% Last update on: 07/12/2015 
% Last use on   : 07/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bDiary = 0;
Diary(mfilename,bDiary);

bSave = 0;
bComparingICRA_speech = 0;
bComparingICRA_piano  = 1;
bICRA_as_paper = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
method = 1; % 1 - script as received
dir = [Get_TUe_paths('ICRA_Tobias')];
addpath([dir delim 'Tools'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bComparingICRA_speech

    fnameS = [Get_TUe_data_paths('db_speechmaterials') 'english' delim 'paul-bagshaw' delim 'wav' delim 'PB-m-all.wav'];
    Nstart=  1;
    Nend  = 5*20000; % 2932000;    

    [signalS fsS] = Wavread(fnameS); 
    
    signalS = signalS(Nstart:Nend);
    lvl = rmsdb(signalS)+100;
    
    noiseS = icra5_noise(signalS,fsS,method);
    % noise2mod = icra5_noise(insig,fs,2);
    % noise2mod = setdbspl(noise2mod,lvl);
    
    % Calculate LTAS
    [ltassdB1S, fHzS] = calcLTAS(signalS,fsS);
    [ltassdB2S      ] = calcLTAS(noiseS,fsS);

    % Calculate MTFs
    [mtf1S,cfModHzS] = calcMTF(signalS,fsS);
    [mtf2S         ] = calcMTF(noiseS,fsS);
    
    %%%
    figure;
    subplot(211);
    semilogx(   fHzS,ltassdB1S,'-b', ... 
                fHzS,ltassdB2S,'-k');
    
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('LTASS [dB]')
    legend({'speech','noise (icra 5)'},'location','southwest')
    axis tight;
    XTick = [20 50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    ylim([-62 -8]);
    
    subplot(212);hold on;
    plot(mtf1S,'.-','marker','s','markerSize',10,'color','b')
    plot(mtf2S,'.-','marker','x','markerSize',10,'color','k')
    set(gca,'xtick',1:numel(cfModHzS),'xticklabel',num2str(cfModHzS'))
    grid on;
    xlabel('Modulation frequency [Hz]')
    ylabel('MTF')
    legend({'speech' 'noise (icra 5)'},'location','northwest')
    
    Save_figure_as(gcf,'ICRA-plot-speech','epsc');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if bComparingICRA_piano
    
    fnameP  = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
    dir_out = [Get_TUe_paths('outputs') 'piano-klanken' delim];
    Mkdir(dir_out);
    
    % fnameP = 'piano-ICRA-tmp.wav';
    [signalP fsP] = Wavread(fnameP); 
    
    lvlP = rmsdb(signalP)+100;
    
    fname_out = sprintf('%spiano_in',dir_out);
    copyfile(fnameP,fname_out);
    
    labels = {'piano-in','icra-5-bands','icra-ERB'};
    
    bCreate = 0;
    %%%
    for i = 2:3 
        
        method    = i; % 2 = 5 banden; 3 = ERB;
        fname_out = sprintf('%spiano-ICRA-method-%.0f',dir_out,method);
        if bCreate
            noiseP    = icra5_noise4piano(signalP,fsP,method);
            Wavwrite(noiseP,fsP,fname_out);
        else
            noiseP  = Wavread(fname_out);
        end
        
        % Calculate LTAS
        [ltassdB1P, fHzP] = calcLTAS(signalP,fsP);
        [ltassdB2P      ] = calcLTAS(noiseP,fsP);

        % Calculate MTFs
        [mtf1P,cfModHzP] = calcMTF(signalP,fsP);
        [mtf2P         ] = calcMTF(noiseP,fsP);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure;
        % subplot(211);
        semilogx(   fHzP,ltassdB1P,'-b', ... 
                    fHzP,ltassdB2P,'-k'); hold on

        grid on;
        xlabel('Frequency [Hz]')
        ylabel('LTASS [dB]')
        legend({labels{1},labels{i}},'location','southwest')
        axis tight;
        XTick = [20 50 100 250 500 1000 2000 4000 8000];
        set(gca,'XTick',XTick);
        set(gca,'XTickLabel',XTick);
        % ylim([-62 -8]);
        % plot([ 800  800], [0 -80],'r--')
        % plot([2400 2400], [0 -80],'r--')
        fname_out = sprintf('%sICRA-plot-piano-LTASS-method-%.0f',dir_out,method);
        Save_figure_as(gcf,fname_out,'epsc');
        
        % subplot(212);hold on;
        figure;
        plot(mtf1P,'.-','marker','s','markerSize',10,'color','b'); hold on
        plot(mtf2P,'.-','marker','x','markerSize',10,'color','k')
        set(gca,'xtick',1:numel(cfModHzP),'xticklabel',num2str(cfModHzP'))
        grid on;
        xlabel('Modulation frequency [Hz]')
        ylabel('MTF')
        legend({labels{1},labels{i}},'location','northwest')

        fname_out = sprintf('%sICRA-plot-piano-MTF-method-%.0f',dir_out,method);
        Save_figure_as(gcf,fname_out,'epsc');
    end
    %%%
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bICRA_as_paper
    
    fnameS_male   = [Get_TUe_data_paths('db_speechmaterials') 'Spanish' delim 'Matrix' delim '00131.wav'];
    fnameS_female = [Get_TUe_data_paths('db_speechmaterials') 'dutch'   delim 'Matrix' delim '00131.wav'];
    
    [signalS fsS] = Wavread(fnameS_male); 
    [signalSf   ] = Wavread(fnameS_female); % assuming the same fsS
    
    Nstart = 1;
    Nend   = min(length(signalS),length(signalSf));
    signalS  = signalS(Nstart:Nend);
    lvl = rmsdb(signalS)+100;
    
    signalSf = signalSf(Nstart:Nend);
    lvlf = rmsdb(signalSf)+100;
    
    method = 1;
    %%%
    opts.Gender = 'male';
    noiseS_n = icra5_noise(signalS         ,fsS,method,'normal',opts); % change normal by raised, loud or shout for other vocal efforts
    label1 = 'male-normal';
    %%%
    opts.Gender = 'female';
    noiseS_r = icra5_noise(signalSf        ,fsS,method,'normal',opts);
    label2 = 'female-normal';
    
    %%%
    opts.Gender = 'male';
    noiseS_l = icra5_noise(signalS+signalSf,fsS,method,'normal',opts);
    label3 = '2-sp normal';
    
    %%%
    % Calculate LTAS
    [ltassdB1S,   fHzS] = calcLTAS( signalS,fsS);
    [ltassdB2S_n      ] = calcLTAS(noiseS_n,fsS);
    [ltassdB2S_r      ] = calcLTAS(noiseS_r,fsS);
    [ltassdB2S_l      ] = calcLTAS(noiseS_l,fsS);

    % Calculate MTFs
    [mtf1S,cfModHzS] = calcMTF(signalS,fsS);
    [mtf2S_n         ] = calcMTF(noiseS_n,fsS); % figure; plot(mtf1S); hold on; plot(mtf2S_n,'r')
    [mtf2S_r         ] = calcMTF(noiseS_r,fsS);
    [mtf2S_l         ] = calcMTF(noiseS_l,fsS);
    
    %%%
    figure;
    subplot(211);
    semilogx(   fHzS,ltassdB2S_n,'-b', ... 
                fHzS,ltassdB2S_r,'-k', ...
                fHzS,ltassdB2S_l,'-r');
    
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('LTASS [dB]')
    legend({label1,label2,label3},'location','southwest')
    axis tight;
    XTick = [20 50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    ylim([-62 0]);
    
    subplot(212);hold on;
    plot(mtf2S_n,'.-','marker','s','markerSize',10,'color','b')
    plot(mtf2S_r,'.-','marker','x','markerSize',10,'color','k')
    plot(mtf2S_l,'.-','marker','o','markerSize',10,'color','r')
    set(gca,'xtick',1:numel(cfModHzS),'xticklabel',num2str(cfModHzS'))
    grid on;
    xlabel('Modulation frequency [Hz]')
    ylabel('MTF')
    legend({label1,label2,label3},'location','northwest')
    
    % Save_figure_as(gcf,'ICRA-plot-speech','epsc');
     
end

rmpath([dir delim 'Tools'])

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

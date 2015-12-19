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
    fnameP = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
    % fnameP = 'piano-ICRA-tmp.wav';
    [signalP fsP] = Wavread(fnameP); 
    
    lvlP = rmsdb(signalP)+100;
    
    method = 3; % ERB
    % noiseP = icra5_noise(signalP,fsP,method);
    noiseP = icra5_noise4piano(signalP,fsP,method);
    % noise2mod = icra5_noise(insig,fs,2);
    % noise2mod = setdbspl(noise2mod,lvl);
    % Wavwrite(noiseP,fsP,'piano-ICRA-tmp');
    
    % Calculate LTAS
    [ltassdB1P, fHzP] = calcLTAS(signalP,fsP);
    [ltassdB2P      ] = calcLTAS(noiseP,fsP);
    
    % Calculate MTFs
    [mtf1P,cfModHzP] = calcMTF(signalP,fsP);
    [mtf2P         ] = calcMTF(noiseP,fsP);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;
    subplot(211);
    semilogx(   fHzP,ltassdB1P,'-b', ... 
                fHzP,ltassdB2P,'-k'); hold on
    
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('LTASS [dB]')
    legend({'speech','noise (icra)'},'location','southwest')
    axis tight;
    XTick = [20 50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    % ylim([-62 -8]);
    plot([ 800  800], [0 -80],'r--')
    plot([2400 2400], [0 -80],'r--')
    
    subplot(212);hold on;
    plot(mtf1P,'.-','marker','s','markerSize',10,'color','b')
    plot(mtf2P,'.-','marker','x','markerSize',10,'color','k')
    set(gca,'xtick',1:numel(cfModHzP),'xticklabel',num2str(cfModHzP'))
    grid on;
    xlabel('Modulation frequency [Hz]')
    ylabel('MTF')
    legend({'piano' 'noise (icra)'},'location','northwest')
        
    Save_figure_as(gcf,'ICRA-plot-piano','epsc');
    
end
rmpath([dir delim 'Tools'])

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

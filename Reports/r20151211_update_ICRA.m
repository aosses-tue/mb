function r20151211_update_ICRA
% function r20151211_update_ICRA
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
bTestICRA       = 0;
bFirstICRApiano = 0; % first version of ICRA piano
bComparingICRA = 1;

if bTestICRA

    file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
    [insig fs] = Wavread(file);
    
    wn = AM_random_noise(20,5000,65,length(insig)/fs,fs);
    outsig = icra5_noise4piano(insig);
    lvl = rmsdb(insig);
    outsig = setdbspl(outsig,lvl+100);
    sound(insig,fs);
    pause(length(insig)/fs*1.2);
    sound(outsig,fs);
    pause(length(insig)/fs*1.2);
    
    outsig2 = icra5_noise4piano(wn);
    outsig2 = setdbspl(outsig2,lvl+100);
    
    sound(outsig2,fs);
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bFirstICRApiano
    
    nRegisters = {'C2','A4','Csh5'};
    
    % Results from loudness balancing:
    Loud_balanc = [ 0   6.6  -2.8; ... % Piano 1 - GRAF 28 (each column a different note)
                    0   0    -4.6; ... % Piano 2 - JBS51-4544
                    0   0     0];      % Piano 3 - NS19
    bFirst_time = 0; % First time means you have not created loudness balanced stimuli yet
    
    for j = 1:length(nRegisters)
        
        dir = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments' delim 'Piano' delim 'Pilot-ICRA' delim 'Stimuli-ICRA-' nRegisters{j} delim]; % 1000-ms-GRAF28-A4_3.wav
        files = Get_filenames(dir,'1000-ms*.wav');
         
        fname = [dir Delete_extension(files{1},'wav')];
        [insigtmp fs] = Wavread([fname '.wav']);
        if bFirst_time
            G = Loud_balanc(1,j);
            Wavwrite(From_dB(G)*insigtmp,fs,[fname '-' num2str(G) '-dB']); % Only first time
        end
        fnoise{1} = [dir files{1+3}];
        insig(:,1) = insigtmp;
         
        fname = [dir Delete_extension(files{2},'wav')]
        [insigtmp fs] = Wavread([fname '.wav']);
        if bFirst_time
            G = Loud_balanc(1,j);
            Wavwrite(From_dB(G)*insigtmp,fs,[fname '-' num2str(G) '-dB']); % Only first time
        end
        fnoise{2} = [dir files{2+3}];
        insig(:,2) = insigtmp;
        
        fname = [dir Delete_extension(files{3},'wav')]
        [insigtmp fs] = Wavread([fname '.wav']);
        if bFirst_time
            G = Loud_balanc(1,j);
            Wavwrite(From_dB(G)*insigtmp,fs,[fname '-' num2str(G) '-dB']); % Only first time
        end
        fnoise{3} = [dir files{3+3}];
        insig(:,3) = insigtmp;

        insig = sum(insig,2);

        outsig = icra5_noise4piano(insig);
        lvl = rmsdb(insig);
        outsig = setdbspl(outsig,lvl+100);

        if bSave
            for i = 1:3
                Wavwrite(outsig,fs,fnoise{i}) % same noise stored three times
            end
        end

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bComparingICRA
    
    % fname = [Get_TUe_paths('ICRA_Tobias') 'EF1_ENG_0001_0.wav'];
    % Nstart=  4000;
    % Nend  = 73000;
    fname = [Get_TUe_data_paths('db_speechmaterials') 'english' delim 'paul-bagshaw-F0mod' delim 'wav' delim 'PB-m-all.wav'];
    Nstart=  1;
    Nend  = 20*20000; % 2932000;    

    % fname = [Get_TUe_data_paths('db_speechmaterials') 'english' delim 'paul-bagshaw-F0mod' delim 'wav' delim 'rl001.wav'];
    % Nstart=  1;
    % Nend  = 30000;    
    
    [insig fs] = Wavread(fname); 
    
    insig = insig(Nstart:Nend);
    lvl = rmsdb(insig)+100;
    
    dir = [Get_TUe_paths('ICRA_Tobias')];
    addpath([dir delim 'Tools'])
    
    if ~isunix % you have to re-compile it for Linux
        noise1 = Demo_Create_Speech_Modulated_Noise(insig,fs);
        noise1 = setdbspl(noise1,lvl);
    end
    
    method = 1; % 1 - script as received
    noise2 = icra5_noise(insig,fs,method);
    noise2mod = icra5_noise(insig,fs,2);
    
    noise2 = setdbspl(noise2,lvl);
    noise2mod = setdbspl(noise2mod,lvl);
    
    speech  = insig;
    noise   = noise2;
    % Calculate LTAS
    [ltassdB1, fHz] = calcLTAS(speech,fs);
    [ltassdB2     ] = calcLTAS(noise,fs);

    % Calculate MTFs
    [mtf1,cfModHz] = calcMTF(speech,fs);
    [mtf2        ] = calcMTF(noise,fs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dur_in = length(insig)/fs;
    
    fname = [Get_TUe_data_paths('db_speechmaterials') 'noise' delim 'icra-misc' delim 'icra3_unmodulated.wav'];
    
    [noise fs] = Wavread(fname); 
    Nstart = fs*2; % starts reading after 2 seconds
    Nend = Nstart + round(dur_in*fs)-1;
    
    noise = noise(Nstart:Nend,1);
    
    % Calculate LTAS
    [ltassdB3, fHz3] = calcLTAS(noise,fs);
    % Calculate MTFs
    [mtf3] = calcMTF(noise,fs);
    
    figure;
    subplot(211);
    semilogx(   fHz,ltassdB1,'-b', ... 
                fHz,ltassdB2,'-k', ...
                fHz3,ltassdB3,'-.r');
    
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('LTASS [dB]')
    legend({'speech','noise (icra 5)','noise (icra 3)'},'location','southwest')
    axis tight;
    XTick = [20 50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    ylim([-62 -8]);
    
    subplot(212);hold on;
    plot(mtf1,'.-','marker','s','markerSize',10,'color','b')
    plot(mtf2,'.-','marker','x','markerSize',10,'color','k')
    plot(mtf3,'.-','marker','o','markerSize',10,'color','r')
    set(gca,'xtick',1:numel(cfModHz),'xticklabel',num2str(cfModHz'))
    grid on;
    xlabel('Modulation frequency [Hz]')
    ylabel('MTF')
    legend({'speech' 'noise (icra 5)' 'noise (icra 3)'},'location','northwest')
    
    Save_figure_as(gcf,'ICRA-plot','epsc');
    
    disp('')
    rmpath([dir delim 'Tools'])
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

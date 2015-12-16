function [noise speech] = Demo_Create_Speech_Modulated_Noise(speech,fsHz)
% function [noise speech] = Demo_Create_Speech_Modulated_Noise(speech,fsHz)

thisfile = mfilename('fullpath');
dir = fileparts([thisfile '.m']);

addpath([dir delim 'Tools'])

if nargin < 1
    fName = 'EF1_ENG_0001_0.wav';
    % Load speech signal
    [speech,fsHz] = Wavread(fName);
end
        
%% MODIFY SIGNAL
% 
% Flatten spectrum while maintaining the modulation content of the signal
noise = processSchroeder(speech);

% Reduce scratchy sound by randomizing the phase
noise = randomizePhase(noise,fsHz);

% Normalize signals
speech = speech / calcRMS(speech);
noise  = noise  / calcRMS(noise);

if nargout == 0
    %% SHOW RESULTS
    % 
    % Calculate LTAS
    [ltassdB1, fHz] = calcLTAS(speech,fsHz);
    [ltassdB2     ] = calcLTAS(noise,fsHz);

    % Calculate MTFs
    [mtf1,cfModHz] = calcMTF(speech,fsHz);
    [mtf2        ] = calcMTF(noise,fsHz);

    figure;
    subplot(211);
    semilogx(fHz,ltassdB1,'-b',fHz,ltassdB2,'-k');
    grid on;
    xlabel('Frequency (Hz)')
    ylabel('LTASS (dB)')
    legend({'speech' 'noise'},'location','southwest')
    axis tight;

    subplot(212);hold on;
    plot(mtf1,'.-','marker','s','markerSize',10,'color','b')
    plot(mtf2,'.-','marker','x','markerSize',10,'color','k')
    set(gca,'xtick',1:numel(cfModHz),'xticklabel',num2str(cfModHz'))
    grid on;
    xlabel('Modulation frequency (Hz)')
    ylabel('MTF')
    legend({'speech' 'noise'},'location','northeast')
end

rmpath([dir delim 'Tools'])

disp('')
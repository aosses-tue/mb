function Generate_Fastl2007_plots
% function Generate_Fastl2007_plots
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 02/09/2014
% Last update on: 02/09/2014 % Update this date manually
% Last use on   : 02/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bFluct = 1; % Chapter 10, fluctuation strength
bRough = 0; % Chapter 11

directory = 'D:\Databases\dir04-Psychoacoustics\Fastl-and-Zwicker-2007\03-Extracted-files\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chapter 10:
if bFluct
    
    % Still offset of +0.1 in values
    filename = [directory 'track_35_t04.wav'];
    x1 = Wavread(filename);
    outs = Do_fluct(x1,44100);
    fmod = 1;
    FS(1) = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod );
    
    filename = [directory 'track_35_t05.wav'];
    x2 = Wavread(filename);
    outs = Do_fluct(x2,44100);
    fmod = 4;
    FS(2) = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod );
    
    filename = [directory 'track_35_t06.wav'];
    x3 = Wavread(filename);
    outs = Do_fluct(x3,44100);
    fmod = 16;
    FS(3) = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod );

    filename = [directory 'track_35_t07.wav'];
    x4 = Wavread(filename);
    outs = Do_fluct(x4,44100);
    fmod = 16;
    FS(4) = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod );
    
    filename = [directory 'track_35_t08.wav'];
    x5 = Wavread(filename);
    outs = Do_fluct(x5,44100);
    fmod = 4;
    FS(5) = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod );
    
    filename = [directory 'track_35_t09.wav'];
    x6 = Wavread(filename);
    outs = Do_fluct(x6,44100);
    fmod = 16;
    FS(6) = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod );
    
    figure;
    plot([1 4 16], FS(1:3), [1 4 16], FS(4:6)), grid on
    legend('AM','FM')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chapter 11:
if bRough
    
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cal file: white noise:
    options.calfile = [directory 'track_03.wav'];
    options.callevel = 60;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure 11.1
    m = [1 0.7 0.4 0.25 0.125 0.1 0];
    
    options.bAnalyser = 12;
    filename = [directory 'track_38_t01.wav'];
    x1 = Wavread(filename);
    m_use(1) = m(1);
    [h ha out1] = PsySoundCL(filename,options);
    
    filename = [directory 'track_38_t03.wav'];
    x2 = Wavread(filename);
    [h ha out2] = PsySoundCL(filename,options);
    m_use(2) = m(3);
    
    filename = [directory 'track_38_t07.wav'];
    x3 = Wavread(filename);
    [h ha out3] = PsySoundCL(filename,options);
    m_use(3) = m(7);
    
    idx = 10;
    figure;
    plot(m_use,[out1.DataRough(idx) out2.DataRough(idx) out3.DataRough(idx)],'o-'), grid on
    
    xlabel('degree of modulation, m')
    ylabel('roughness [asper]')
    title(sprintf('f_c = 1kHz, fmod = 70 Hz, SPL = %.2f [dB]',rmsdb(x3)+90))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

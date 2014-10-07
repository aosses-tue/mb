% Meas 2013 01 28:
%
% Dependencies:
%
%   name2figname
%   Electrodogram_Sync
%   Subplot_sequence
%   in22channel
%   AverageAmplitudePerChannel
%   FigureAverageChannel

clear, clc, close all

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('Begin script "meas_20130128.m"')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference electrodogram (Choose one and then truncate it!)

Fs = 16000;

FormatElectrodogram.loop = 'n'; % Looped recording
FormatElectrodogram.CompChannel = 18;
FormatElectrodogram.Length_analysis = 1; % Time in seconds
FormatElectrodogram.Path = '/home/alejandro/Documenten/Meas/20130128_Cal_SP12';

k = 55:65;
% k = [k 70:6:82];

for i = k
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ACE Results:
%     name1 = 'cal_55dB_20130125.out';
    name1 = ['cal_' num2str(i) 'dB_20130128.out'];

    % F0 Results:
    name2 = ['cal_' num2str(i) 'dB_20130128.out'];

    [p{i}.pACE, p{i}.pF0m] = Electrodogram_Sync(name1,name2, FormatElectrodogram);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    y  = in22channel(p{i}.pACE.electrodes, p{i}.pACE.current_levels);
    
    p{i}.pACE.Energy = AverageAmplitudePerChannel(y,FormatElectrodogram.Length_analysis,Fs,9);
    
    y  = in22channel(p{i}.pF0m.electrodes, p{i}.pF0m.current_levels);

    p{i}.pF0m.Energy = AverageAmplitudePerChannel(y,FormatElectrodogram.Length_analysis,Fs,9);
    
    MatrixACE(i,1:22)   = p{i}.pACE.Energy;
    
    MatrixMeanACE(i)    = ( p{i}.pACE.Energy( 1,7 ) ); 
    
    MatrixF0m(i,1:22)   = p{i}.pF0m.Energy;
    
    % name1Fig    = name2figname(name1);
    % name2Fig    = name2figname(name2);

    % Energ = AverageAmplitudePerChannel(y,6,Fs,9);
    % Energ0 = AverageAmplitudePerChannel(y0,6,Fs,9);
    % FigureAverageChannel(transpose([Energ; Energ0])/250,name1Fig,name2Fig), title('Amplitudes Normalised to 250')
end

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('End script "meas_20130128.m"')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
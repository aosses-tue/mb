function PsySoundCL(filename)
% function PsySoundCL(filename)
%   
%   Executes PsySound3 from command line
% 
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 22/7/2014
% Last update on: 29/7/2014 % Update this date manually
% Last use on   : 29/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1:

% Daniel and Weber Roughness model
close all
if nargin == 0
    % filename = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\modus-1_v2-2filt-fc-1000-Hz.wav';
    filename = 'D:\Documenten-TUe\10-Referenties\02-Mijn-boeken\Zwicker-Psychoacoustics\Sound\track_42.wav';
end

[x fs] = Wavread(filename);
% rmsdb(x)
dB_value = dbspl(x);
% 0 dBFS = 100 dB

fh = readData(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Calibation
% % To calibrate your files by using a file ‘CalFile.wav’ that represents 91 dB:
% fhs = calibrate(fhs, ’WithFiles’, ’CalFile.wav’, 91)
fh = calibrate(fh, 'WithFiles', filename, dB_value); % calibrated with itself
% fh.calCoeff = 1; % To avoid calibration set this value to 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Analysers

% obj = SLM(fh); 
obj = RoughnessDW(fh);
obj = process(obj,fh,[]);

tmpObj  = get(obj,'output');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: PostProcessing

% Average Roughness -> Graph (Visualisation - Single Axis, line)
DataSpecRoughness = get(tmpObj{1,2},'Data');
z       = get(tmpObj{1,2},'Freq');
figure
plot(z, mean(DataSpecRoughness) );
xlabel('Critical band rate (Bark)')
ylabel('Specific Roughness (Aspers/Bark)')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roughness -> Graph (Visualisation - Single Axis, plot)
t       = get(tmpObj{1,1},'Time');
DataRoughness = get(tmpObj{1,1},'Data');
figure;
plot(t, DataRoughness);
xlabel('Time (seconds)')
ylabel('Roughness (aspers)')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % To estimate the time the analysis will take
% runanalysis(fhs, ’FFT’, ’Hilbert’, ’SLM’, ’estimate’)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
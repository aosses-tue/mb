function [h ha] = PsySoundCL(filename,option)
% function [h ha] = PsySoundCL(filename,option)
% 
% 1. Description:
%       Executes PsySound3 from command line
%       Specify:
%           option.nAnalyser    - 12 = DLM Chalupper's model
%                               - 15 = Daniel's Roughness
%           option.calfile
%           option.callevel
%   
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       PsySoundCL;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 22/07/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
ha = [];

% Step 1:

% Daniel and Weber Roughness model
% close all
if nargin < 2
    option = [];
end

option = Ensure_field(option,'nAnalyser',12);

nAnalyser = option.nAnalyser;

if nargin == 0
    % filename = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\modus-1_v2-2filt-fc-1000-Hz.wav';
    % filename = 'D:\Documenten-TUe\10-Referenties\02-Mijn-boeken\Zwicker-Psychoacoustics\Sound\track_42.wav'; % 78.7013 dB
    filename = 'D:\Documenten-TUe\10-Referenties\02-Mijn-boeken\Fastl2007-psychoacoustics\Extracted\ref_sharpness.wav';
end

if ~isfield(option,'title')
    option = Ensure_field(option,'title',name2figname(filename));
end
option = Ensure_field(option,'colorbar_scale','Gray');

if nargin ~= 1 & isnumeric(filename)
    x = filename;
    fs = option.fs;
else
    [x fs] = Wavread(filename);
end

if ~isfield(option,'calfile')
    disp([mfilename '.m: calibration file not specified. 0 dBFS = 100 dB SPL is going to be used...'])
    option.calfile = filename;
    option.callevel = dbspl(x);
else
    if ~isfield(option,'callevel')
        option.callevel = input([mfilename '.m - Specify the level of the calibration tone [dB SPL] : ']);
    end
end

fh = readData(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Calibation
% % To calibrate your files by using a file 'CalFile.wav' that represents 91 dB:
% fhs = calibrate(fhs, 'WithFiles', 'CalFile.wav', 91)
fh = calibrate(fh, 'WithFiles', option.calfile, option.callevel); 
% fh.calCoeff = 1; % To avoid calibration set this value to 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Analysers

switch nAnalyser
    case 8
        obj = SLM(fh); 
    case 12
        obj = LoudnessCF(fh);
    case 15
        obj = RoughnessDW(fh);
end

obj = process(obj,fh,[]);
tmpObj  = get(obj,'output');

t       = get(tmpObj{1,1},'Time');
z       = get(tmpObj{1,2},'Freq'); % 1:24

switch nAnalyser 
    case 12
        
        zspec = get(tmpObj{1,5},'Freq');
        
        for i = 1:length(tmpObj)
            disp( get(tmpObj{1,i},'Name') );
        end
        
        % Average main loudness:
        nParam = 4;
        Data = get(tmpObj{1,nParam},'Data');
        
        figure;
        plot(z,Data);
        xlabel('Critical band rate (Bark)')
        ylabel('Loudness (Sones/Bark)');
        title(sprintf('Average Main Loudness - %s', option.title));
        grid on;
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        % Average specific loudness
        nParam = 5;
        DataAvSpecLoud = get(tmpObj{1,nParam},'Data');
        
        figure;
        plot(zspec,DataAvSpecLoud);
        xlabel('Critical band rate (Bark)')
        ylabel('Loudness (Sones/Bark)');
        title(sprintf('Average Specific Loudness - %s', option.title));
        grid on;
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        % Loudness
        nParam = 1;
        DataLoud = get(tmpObj{1,nParam},'Data');
        
        figure;
        plot(t,DataLoud)
        xlabel('Time (Seconds)')
        ylabel('Loudness (Sones)');
        title(sprintf('Loudness - %s', option.title));
        grid on;
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        % Main loudness
        nParam = 2;
        DataMainLoud = get(tmpObj{1,nParam},'Data');
        
        figure;
        imagesc(t,z,DataMainLoud');
        set(gca,'YDir','Normal');
        colormap(option.colorbar_scale);
        hcb = colorbar; 
        ylabel('Critical band rate (Bark)')
        xlabel('Time (seconds)')
        set(get(hcb,'YLabel'),'String','Loudness (Sones/Bark)')
        title(sprintf('Main loudness - %s', option.title));
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        % Specific loudness
        nParam = 3;
        DataSpecLoud = get(tmpObj{1,nParam},'Data');
        
        figure;
        imagesc(t,zspec,DataSpecLoud');
        set(gca,'YDir','Normal');
        colormap(option.colorbar_scale);
        hcb = colorbar; 
        ylabel('Critical band rate (Bark)')
        xlabel('Time (seconds)')
        set(get(hcb,'YLabel'),'String','Loudness (Sones/Bark)')
        title(sprintf('Specific loudness - %s', option.title));
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        nParam = 6;
        DataSharp = get(tmpObj{1,nParam},'Data');
        
        figure;
        plot(t,DataSharp)
        xlabel('Time (Seconds)')
        ylabel('Sharpness (Acums)');
        title(sprintf('Sharpness - %s', option.title));
        grid on;
        h(end+1) = gcf;
        ha(end+1) = gca;
        
    case 15
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: PostProcessing

        % Average Roughness -> Graph (Visualisation - Single Axis, line)
        DataSpecRoughness = get(tmpObj{1,2},'Data');
        
        figure
        plot(z, mean(DataSpecRoughness) );
        xlabel('Critical band rate (Bark)')
        ylabel('Specific Roughness (Aspers/Bark)')
        title(sprintf('Average Roughness - %s', option.title));
        grid on
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Roughness -> Graph (Visualisation - Single Axis, plot)
        
        DataRoughness = get(tmpObj{1,1},'Data');
        figure;
        plot(t, DataRoughness);
        xlabel('Time (seconds)')
        ylabel('Roughness (aspers)')
        title(sprintf('Roughness - %s', option.title));
        grid on
        
        h(end+1) = gcf;
        ha(end+1) = gca;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% % To estimate the time the analysis will take
% runanalysis(fhs, �FFT�, �Hilbert�, �SLM�, �estimate�)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
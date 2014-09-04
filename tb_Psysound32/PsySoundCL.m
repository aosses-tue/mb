function [h ha output] = PsySoundCL(filename,option)
% function [h ha output] = PsySoundCL(filename,option)
% 
% 1. Description:
%       Executes PsySound3 from command line
%       Analysers tested:
%            8 - SLM(fh); 
%           10 - ThirdOctaveBand(fh);   % NOT VALIDATED YET
%           11 - CPBFFT(fh);            % NOT VALIDATED YET
%           12 - LoudnessCF(fh);
%           15 - RoughnessDW(fh);
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       PsySoundCL;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 22/07/2014
% Last update on: 29/08/2014 % Update this date manually
% Last use on   : 29/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
ha = [];

% Step 1:

% Daniel and Weber Roughness model
% close all
if nargin < 2
    option = [];
end

if nargin == 0
    option.nAnalyser = input('Choose analyser to be used: \n - type 12 for DLM model or \n - type 15 for Roughness model :');
    [f1 f2] = uigetfile(Get_TUe_paths('outputs'));
    filename = [f2 f1];
end

option = Ensure_field(option,'nAnalyser'  ,15);
option = Ensure_field(option,'bCosineRamp',0 );

if option.bCosineRamp
    option = Ensure_field(option,'bCosineRampOnset',25);
end

nAnalyser = option.nAnalyser;
% nAnalyser = 12; % Dynamic loudness
% nAnalyser = 15; % Roughness

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
    disp([mfilename '.m: calibration file not specified. Calibration respect to itself, 0 dBFS = 100 dB SPL is going to be used...'])
    
    bCal = input('Choose calibration method: (1-AMT for 0 dBFS = 100 dB / 2-Zwicker for 0 dBFS = 90 dB): ');
    
    switch bCal
        case 1
            option.calfile = filename;
            option.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav']; %white noise, adjusted to AMT convention
            option.callevel = 70;
        case 2
            option.calfile = [Get_TUe_paths('db_fastl2007') 'track_03.wav']; %white noise
            option.callevel = 60;
    end
else
    if ~isfield(option,'callevel')
        option.callevel = input([mfilename '.m - Specify the level of the calibration tone [dB SPL] : ']);
    end
end

if ~option.bCosineRamp 
    fh = readData(filename); % then no ramp is applied
else
    ramp2apply = cos_ramp(length(x),fs,option.bCosineRampOnset);
    x = transpose(ramp2apply).*x;
    
    [xx filenametmp] = fileparts(filename); 
    dirtmp = [Get_TUe_paths('outputs') delim 'tmp-PsySound' delim];
    filenametmp = [dirtmp filenametmp '.wav'];
    Mkdir(dirtmp);
    Wavwrite(x,fs,filenametmp);
    fh = readData(filenametmp);
    fprintf('file calibrated to %.2f [dB SPL] if Zwicker''s tone were used\n',rmsdb(x)+90);
    fprintf('file calibrated to %.2f [dB SPL] if Zwicker''s tone were used if AMT values are used\n',rmsdb(x)+100);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Calibation
disp(['Calibration file      : ' option.calfile])
disp(['Calibration level [dB]: ' num2str(option.callevel)])
fh = calibrate(fh, 'WithFiles', option.calfile, option.callevel); 
% fh.calCoeff = 1; % To avoid calibration set this value to 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Analysers

switch nAnalyser
    case 8
        obj = SLM(fh); 
    case 10
        obj = ThirdOctaveBand(fh);
    case 11
        obj = CPBFFT(fh);
    case 12
        obj = LoudnessCF(fh);
    case 15
        obj = RoughnessDW(fh);
end

obj = process(obj,fh,[]);
tmpObj  = get(obj,'output');

t       = get(tmpObj{1,1},'Time');

if nAnalyser == 12 | nAnalyser == 15
    
    z   = get(tmpObj{1,2},'Freq'); % 1:24
    
elseif nAnalyser == 10 | nAnalyser == 11
    
    f_cell = get(tmpObj{1,2},'Freq');
    f = zeros(size(f_cell));
    for i=1:length(f_cell)
        f(i) = str2num(f_cell{i});
    end
    
    f_cell = get(tmpObj{1,3},'Freq');
    f_octv = zeros(size(f_cell));
    for i=1:length(f_cell)
        f_octv(i) = str2num(f_cell{i});
    end
end

if isfield(option,'trange')
    % for i = 1:size(option.trange,1)
    %     tmp = Convert2nan_if_outofrange(t,option.trange);
    % end
    t = option.trange;
end
    
output.t = t;
output.z = z;

switch nAnalyser 
    
    case 10
        % Out1
        % Out2
        % Out3: Loudness with percentiles
        
        
    case 11
        
        % One-third Octave Band Spectrogram
        % Size = Number of frames x 33 bands(for 1/3 OB)
        DataSpecOneThird = get(tmpObj{1,1},'Data');
        
        % One-third Octave Band Spectrogram, global values
        DataSpecOneThirdAvg = get(tmpObj{1,2},'Data');
        
        DataSpecOctv    = get(tmpObj{1,3},'Data');
        DataSpecOctvAvg = get(tmpObj{1,4},'Data');
        
        figure;
        semilogx(f,DataSpecOneThirdAvg,'o-');
        ylabel('Magnitude (dB)')
        xlabel('Frequency (Hz)');
        title(sprintf('One-third octave band spectrum - %s', option.title));
        grid on;
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        output.f        = f;
        output.f_octv   = f_octv;
        output.DataSpecOneThird     = DataSpecOneThird;             
        output.DataSpecOneThirdAvg  = DataSpecOneThirdAvg; 
        output.DataSpecOctv         = DataSpecOctv;             
        output.DataSpecOctvAvg      = DataSpecOctvAvg;  
        
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
        
        output.zspec = zspec;
        
        output.DataLoud = DataLoud;             % Param 1: Loudness
                                                % Param 2: Main loudness - 3D
                                                % Param 3: Specific loudness - 3D
                                                % Param 4: Average Main loudness % not interesting by now
        output.DataAvSpecLoud = DataAvSpecLoud; % Param 5: Average specific loudness
        output.DataSharp      = DataSharp;      % Param 6: Sharpness
        
    case 15
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: PostProcessing

        % Average Roughness -> Graph (Visualisation - Single Axis, line)
        DataSpecRough = get(tmpObj{1,2},'Data');
        
        figure
        plot(z, mean(DataSpecRough) );
        xlabel('Critical band rate (Bark)')
        ylabel('Specific Roughness (Aspers/Bark)')
        title(sprintf('Average Roughness - %s', option.title));
        grid on
        h(end+1) = gcf;
        ha(end+1) = gca;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Roughness -> Graph (Visualisation - Single Axis, plot)
        
        DataRough = get(tmpObj{1,1},'Data');
        figure;
        plot(t, DataRough);
        xlabel('Time (seconds)')
        ylabel('Roughness (aspers)')
        title(sprintf('Roughness - %s', option.title));
        grid on
        
        h(end+1) = gcf;
        ha(end+1) = gca;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        output.DataRough = DataRough;
        output.DataSpecRough = DataSpecRough;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
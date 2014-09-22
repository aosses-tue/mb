function [output h ha] = PsySoundCL(filename,option)
% function [output h ha] = PsySoundCL(filename,option)
% 
% PAS OP:
%       before 17/09: [h ha output] = PsySoundCL(filename,option)
%       after 17/09: [output h ha] = PsySoundCL(filename,option)
% 
% 1. Description:
%       Executes PsySound3 from command line
% 
%       filename (Include extension) - file to be analysed
%       option.CalMethod - file calibration according to known references
%       option.nAnalyser - analyser to be applied. Analysers tested:
%            8 - SLM(fh); 
%           10 - ThirdOctaveBand(fh);
%           11 - CPBFFT(fh);
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
% Last update on: 17/09/2014 % Update this date manually
% Last use on   : 17/09/2014 % Update this date manually
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
    str_analysers = ['\n - type 10 for one-third OB or ' ...
                     '\n - type 12 for DLM model or ' ...
                     '\n - type 15 for Roughness model'];
                
    option.nAnalyser = input(['Choose analyser to be used: ' str_analysers ' :']);
    [f1 f2] = uigetfile(Get_TUe_paths('outputs'));
    filename = [f2 f1];
end

option = Ensure_field(option,'bPlot'      , 1)
option = Ensure_field(option,'nAnalyser'  ,15);
option = Ensure_field(option,'bCosineRamp', 0);

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
    
    if ~isfield(option,'CalMethod')
        bCal = input('Choose calibration method: (1-AMT for 0 dBFS = 100 dB / 2-Zwicker for 0 dBFS = 90 dB): ');
    else
        bCal = option.CalMethod;
    end
        
    switch bCal
        case 1
            option.calfile = filename;
            option.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav']; %white noise, adjusted to AMT convention
            option.callevel = 70;
        case 2
            option.calfile = [Get_TUe_paths('db_fastl2007') 'track_03.wav']; %white noise
            option.callevel = 60;
        case 3
            tmp = Get_TUe_subpaths('db_speechmaterials');
            option.calfile = [tmp.allfiles_PB 'whitenoise-f.wav']; % white noise
            option.callevel = 65;
        case 4
            tmp = Get_TUe_subpaths('db_speechmaterials');
            option.calfile = [tmp.allfiles_PB 'whitenoise-m.wav']; % white noise
            option.callevel = 65;
        case 5
            tmp = Get_TUe_subpaths('db_speechmaterials');
            option.calfile = [tmp.allfiles_LISTf 'wivineruis.wav']; % SSN
            option.callevel = 65;
    end
else
    if ~isfield(option,'callevel')
        option.callevel = input([mfilename '.m - Specify the level of the calibration tone [dB SPL] : ']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO: work on something similar to readData but using data 'not stored'

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
        option = Ensure_field(option,'Author','DH');
        
        obj = RoughnessDW(fh);
        
        if strcmp(option.Author,'DH') % Dik Hermes
            
            st.type = 'samples'; 
            st.size = 4096; 
            obj = set(obj,'overlap',st);
            
        end
        fprintf('Roghness model by %s is being used. Default is DW\n',option.Author);
        pause(2)
    case 20
        
        obj = FluctuationStrength(fh);
        st.type = 'samples'; 
        st.size = 4096; 
        obj = set(obj,'overlap',st);
        
end

obj = process(obj,fh,[]);
tmpObj  = get(obj,'output');

t       = get(tmpObj{1,1},'Time');

if nAnalyser == 12 | nAnalyser == 15 | nAnalyser == 20
    
    z   = get(tmpObj{1,2},'Freq'); % 1:24
    
elseif nAnalyser == 10
    
    f_cell = get(tmpObj{1,2},'Freq');
    f = zeros(size(f_cell));
    for i=1:length(f_cell)
        f(i) = str2num(f_cell{i});
    end
    
elseif nAnalyser == 11 
    
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

if exist('z','var')
    output.z = z;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stats Nparams] = Get_stats_PsySound(tmpObj);

switch nAnalyser 
    
    case 8
        % Out1: SPL A-weighted Slow
        % Out2: SPL A-weighted Slow
        % Out3: SPL Z-Unweighted Fast
        % Out4: SPL Z-Unweighted Slow
        % Nparams = 4;
        
        nParam = 1;
        Data_dBAF = get(tmpObj{1,nParam},'Data');
        nParam = 2;
        Data_dBAS = get(tmpObj{1,nParam},'Data');
        nParam = 3;
        Data_dBZF = get(tmpObj{1,nParam},'Data');
        nParam = 4;
        Data_dBZS = get(tmpObj{1,nParam},'Data');
        output.description = {'dB(A) Fast','dB(A) Slow','dB(Z) Fast','dB(Z) Slow'};
        
        output.Data_dBAF = Data_dBAF;
        output.Data_dBAS = Data_dBAS;
        output.Data_dBZF = Data_dBZF;
        output.Data_dBZS = Data_dBZS;
        
        
    case 10
        % Out1: 
        % Out2
        % Out3: Loudness with percentiles
        nParam = 1;
        DataSpecOneThird = get(tmpObj{1,nParam},'Data');        
        
        nParam = 2;
        DataSpecOneThirdAvg = get(tmpObj{1,nParam},'Data');
        txt_title = sprintf('%s - %s',get(tmpObj{1,nParam},'Name'),option.title);
        
        if option.bPlot
            figure;
            semilogx(f,DataSpecOneThirdAvg,'o-');
            ylabel('Magnitude (dB)')
            xlabel('Frequency (Hz)');
            title(txt_title); % title(sprintf('One-third octave band spectrum - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        nParam = 3;
        zspec = transpose( get(tmpObj{1,nParam},'Freq') ); 
        Data  = get(tmpObj{1,nParam},'Data');
        txt_title = sprintf('%s - %s',get(tmpObj{1,nParam},'Name'),option.title);
        
        if option.bPlot
            figure;
            plot(zspec,Data);
            xlabel('Critical band rate (Bark)')
            ylabel('Loudness (Sones/Bark)');
            title(txt_title); % title(sprintf('Specific Loudness (ISO532B) - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        output.f = f;
        output.DataSpecOneThird     = DataSpecOneThird;
        output.DataSpecOneThirdAvg  = DataSpecOneThirdAvg;
        
        output.zspec = zspec;
        output.DataLoud = Data;
        stats.loud_tot = sum(Data)*0.1;
        
        % output.percentilesLoud
        
    case 11
        
        % One-third Octave Band Spectrogram
        % Size = Number of frames x 33 bands(for 1/3 OB)
        nParam = 1;
        DataSpecOneThird = get(tmpObj{1,nParam},'Data');
        
        % One-third Octave Band Spectrogram, global values
        nParam = 2;
        DataSpecOneThirdAvg = get(tmpObj{1,nParam},'Data');
        
        nParam = 3;
        DataSpecOctv    = get(tmpObj{1,nParam},'Data');
        
        nParam = 4;
        DataSpecOctvAvg = get(tmpObj{1,nParam},'Data');
        
        % output.description = {'','','',''};
        
        if option.bPlot
            figure;
            semilogx(f,DataSpecOneThirdAvg,'o-');
            ylabel('Magnitude (dB)')
            xlabel('Frequency (Hz)');
            title(sprintf('One-third octave band spectrum - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
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
        
        if option.bPlot
            figure;
            plot(z,Data);
            xlabel('Critical band rate (Bark)')
            ylabel('Loudness (Sones/Bark)');
            title(sprintf('Average Main Loudness - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        % Average specific loudness
        nParam = 5;
        DataAvSpecLoud = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
            figure;
            plot(zspec,DataAvSpecLoud);
            xlabel('Critical band rate (Bark)')
            ylabel('Loudness (Sones/Bark)');
            title(sprintf('Average Specific Loudness - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        % Loudness
        nParam = 1;
        DataLoud = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
            figure;
            plot(t,DataLoud)
            xlabel('Time (Seconds)')
            ylabel('Loudness (Sones)');
            title(sprintf('Loudness - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        % Main loudness
        nParam = 2;
        DataMainLoud = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
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
        end
        
        % Specific loudness
        nParam = 3;
        DataSpecLoud = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
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
        end
        
        nParam = 6;
        DataSharp = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
            figure;
            plot(t,DataSharp)
            xlabel('Time (Seconds)')
            ylabel('Sharpness (Acums)');
            title(sprintf('Sharpness - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
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

        DataRough       = get(tmpObj{1,1},'Data');
        DataSpecRough   = get(tmpObj{1,2},'Data');
        
        if option.bPlot
            figure
            plot(z, mean(DataSpecRough) );
            xlabel('Critical band rate (Bark)')
            ylabel('Specific Roughness (Aspers/Bark)')
            title(sprintf('Average Roughness - %s', option.title));
            grid on
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Roughness -> Graph (Visualisation - Single Axis, plot)
        
        if option.bPlot
            figure;
            plot(t, DataRough);
            xlabel('Time (seconds)')
            ylabel('Roughness (aspers)')
            title(sprintf('Roughness - %s', option.title));
            grid on

            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        output.DataRough = DataRough;
        output.DataSpecRough = DataSpecRough;
        stats.rough_tot = 0.5*sum(  mean( DataRough )  );
        
	case 20
        
        DataFluct       = get(tmpObj{1,1},'Data');
        DataSpecFluct   = get(tmpObj{1,2},'Data');
        
        if option.bPlot
            figure
            plot(z, mean(DataSpecFluct) );
            xlabel('Critical band rate (Bark)')
            ylabel('Specific Fluctuation Strength (Vacils/Bark)')
            title(sprintf('Average Fluctuation Strength - %s', option.title));
            grid on
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Roughness -> Graph (Visualisation - Single Axis, plot)
       
        if option.bPlot
            figure;
            plot(t, DataFluct);
            xlabel('Time (seconds)')
            ylabel('Fluctuation Strength (vacils)')
            title(sprintf('Fluctuation Strength - %s', option.title));
            grid on
        
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        output.DataFluct = DataFluct;
        output.DataSpecFluct = DataSpecFluct;
        stats.fluct_tot = 0.5*sum(  mean( DataFluct )  );
        disp('')
end

output.nAnalyser = nAnalyser;
output.stats     = stats;

if nargout == 0
    
    
    for i = 1:Nparams
        fprintf('Analyser Number: %.0f\n',nAnalyser);
        fprintf('Param %.0f out of %.0f: %s\n',i,Nparams,get(tmpObj{1,i},'Name') );
        try
            fprintf('\t [average] = [%.3f]\n'     , stats.mean(i));
            fprintf('\t [min max] = [%.3f %.3f]\n', stats.min(i),stats.max(i));
        catch
            fprintf('\t It seems that this parameter corresponds to a Spectrum object\n');
        end
    end
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
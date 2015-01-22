function [output h ha] = PsySoundCL(filename,option,st)
% function [output h ha] = PsySoundCL(filename,option,st)
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
%           option.nAnalyser = 1;  %  
%           option.nAnalyser = 8;  % SLM(fh); 
%           option.nAnalyser = 10; % ThirdOctaveBand(fh);
%           option.nAnalyser = 11; % CPBFFT(fh);
%           option.nAnalyser = 12; % LoudnessCF(fh);
%           option.nAnalyser = 15; % RoughnessDW(fh);
%           option.nAnalyser = 20; % FluctuationStrength(fh);
%       st - to pass config parameters to object
% 
% 2. Stand-alone example:
%       PsySoundCL;
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% 4. TO DO:
%       work on something similar to readData but using data 'not stored' (L114)

% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 22/07/2014
% Last update on: 20/01/2015 % Update this date manually
% Last use on   : 20/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
ha = [];

% close all
if nargin < 2
    option = [];
end

if nargin == 0
    [f1 f2] = uigetfile(Get_TUe_paths('outputs'));
    filename = [f2 f1];
end

if nargin < 2
    str_analysers = ['\n - type 1, 8' ...
                     '\n - type 10 for one-third OB or ' ...
                     '\n - type 12 for DLM model or ' ...
                     '\n - type 15 for Roughness model'];
                
    option.nAnalyser = input(['Choose analyser to be used: ' str_analysers ' :']);
end

option = Ensure_field(option,'bPlot'      , 1);
option = Ensure_field(option,'nAnalyser'  ,1);

nAnalyser = option.nAnalyser;

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

%% Calibration, reference file:
%       Selects reference file to adjust proper SPL level for input signals
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
        % case 99 % custom
        %     option.calfile = filename;
        %     option.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav']; %white noise, adjusted to AMT convention
        %     option = Ensure_field(option,'callevel',70);
    end
else
    if ~isfield(option,'callevel')
        option.callevel = input([mfilename '.m - Specify the level of the calibration tone [dB SPL] : ']);
    end
end

%% Extracts excerpt:
if isfield(option,'tanalysis') 
    
    option = Ensure_field(option,'bGenerateExcerpt',1);
    try
        filename_excerpt = [Delete_extension(filename,'wav') '-e.wav']; 
        
        if option.bGenerateExcerpt == 1
            [xx ffs] = Wavread(filename);
            tt = (1:length(xx))/ffs;
            idx = find(tt>=option.tanalysis(1) & tt<=option.tanalysis(2));
    
            Wavwrite( xx(idx),ffs,filename_excerpt );
            filename = filename_excerpt;
        else
            [xx ffs] = Wavread(filename);
            
            Wavwrite( xx(idx),ffs,filename_excerpt );
            filename = filename_excerpt;
        end
    catch
        option.bGenerateExcerpt = 0;
    end
    
    if option.bGenerateExcerpt == 0
        warning('Excerpt not extracted. Analysis will be done using complete audio file');
    end
    
    fh = readData(filename); % then no ramp is applied
    
else
    fh = readData(filename);
end    

%% Step 2: Calibation: applying adjustments
disp(['Calibration file      : ' option.calfile])
disp(['Calibration level [dB]: ' num2str(option.callevel)])
fh = calibrate(fh, 'WithFiles', option.calfile, option.callevel); 
% fh.calCoeff = 1; % To avoid calibration set this value to 1

%% Step 3: Analysers

switch nAnalyser
    case 1
        obj = FFT(fh);
        
        obj = set(obj,'windowLength',2*8192);
        st = [];
        st.type = 'percent';
        st.size = 75; % default = 75
        obj = set(obj,'overlap',st);

    case 8
        obj = SLM(fh); 
    case 10
        obj = ThirdOctaveBand(fh);
    case 11
        obj = CPBFFT(fh);
        
        obj = set(obj,'windowLength',8192);
        st = [];
        st.type = 'percent';
        st.size = 75; % default = 75
        obj = set(obj,'overlap',st);
        
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
        st.type = 'samples'; % commented temporarily on 06/11/2014
        
        tmp = get(obj,'samples');

        % st.type = 'samples'; 
        % st.size = 4096; 
        % obj = set(obj,'overlap',st);
        
end

obj = process(obj,fh,[]);
tmpObj  = get(obj,'output');

t       = get(tmpObj{1,1},'Time');

if nAnalyser == 1
    
    f   = get(tmpObj{1,2},'Freq'); % 1:fs/2
    
elseif nAnalyser == 12 | nAnalyser == 15 | nAnalyser == 20
    
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
    t = option.trange;
end
    
output.t = t;

if exist('z','var')
    output.z = z;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stats Nparams] = Get_stats_PsySound(tmpObj);

switch nAnalyser 
    
    case 1
        
        nParam = 1;
        Data1   = get(tmpObj{1,nParam},'Data');
        output.Data1 = Data1;
        output.name{nParam} = get(tmpObj{1,nParam},'Name');
        output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');
        
        nParam = 2;
        Data2   = get(tmpObj{1,nParam},'Data');
        % Data2name = '';
        output.Data2 = Data2;
        output.name{nParam} = get(tmpObj{1,nParam},'Name');
        output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');
        output.f = f;
        
        if option.bPlot
            figure;
            semilogx(f,Data2,'-');
            xlabel('Frequency [Hz]');
            ylabel('Amplitude [dB]');
            grid on
        end 
        
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
        DataAvMainLoud  = get(tmpObj{1,nParam},'Data');
        txt_title = sprintf('%s - %s',get(tmpObj{1,nParam},'Name'),option.title);
        
        if option.bPlot
            figure;
            plot(zspec,DataAvMainLoud);
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
        output.DataLoud = DataAvMainLoud;
        stats.loud_tot = sum(DataAvMainLoud)*0.1;
        
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
        DataAvMainLoud = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
            figure;
            plot(z,DataAvMainLoud);
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
        Data1 = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
            figure;
            plot(t,Data1)
            xlabel('Time (Seconds)')
            ylabel('Loudness (Sones)');
            title(sprintf('Loudness - %s', option.title));
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
        end
        
        % Main loudness
        nParam = 2;
        Data2 = get(tmpObj{1,nParam},'Data');
        
        if option.bPlot
            figure;
            imagesc(t,z,Data2');
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
        output.DataLoud       = Data1;          % Param 1: Loudness
        output.DataMainLoud   = Data2;          % Param 2: Main loudness - 3D
        output.DataSpecLoud   = DataSpecLoud;   % Param 3: Specific loudness - 3D
        output.DataAvMainLoud = DataAvMainLoud; % Param 4: Average Main loudness % not interesting by now
        output.DataAvSpecLoud = DataAvSpecLoud; % Param 5: Average specific loudness
        output.DataSharp      = DataSharp;      % Param 6: Sharpness
        
    case 15
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: PostProcessing
        
        nParam  = 1;
        Data1   = get(tmpObj{1,nParam},'Data');
        output.name{nParam}     = get(tmpObj{1,nParam},'Name');
        output.param{nParam}    = strrep( lower( output.name{nParam} ),' ','-');
        
        nParam  = 2;
        Data2   = get(tmpObj{1,nParam},'Data');
        output.name{nParam}     = get(tmpObj{1,nParam},'Name');
        output.param{nParam}    = strrep( lower( output.name{nParam} ),' ','-');
        
        output.Data1    = Data1;
        output.Data2    = Data2;
        
        DataRough       = Data1;
        DataSpecRough   = Data2;
        
        stats.rough_tot = mean( DataRough );
        
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
        stats.fluct_tot = mean( DataFluct );
        % stats.fluct_tot = mean( 0.25*sum( DataSpecFluct' ) );
        disp('')
end

output.nAnalyser = nAnalyser;
output.stats     = stats;

if option.bGenerateExcerpt
    
    disp('...deleting temporal audio file');
    delete( filename_excerpt );
    
else
    
    disp('...deleting temporal audio file, with zero padding');
    delete( filename_excerpt );
    
end

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
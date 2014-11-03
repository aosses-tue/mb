function [h ha stats] = PsySoundCL_Figures_3(param,res1, res2, res3, option, stPlot)
% function [h ha stats] = PsySoundCL_Figures_3(param,res1, res2, res3, option, stPlot)
% 
% 1. Description:
%       To plot PsySound results for 2 data series
%       'param' can be:
% 
%       nAnalyser                           excerpt     stats
%       1           - 'FFT'                 NO          NO
%       10          - 'one-third-OB'        NO          NO
%       10          - 'specific-loudness'   NO          NO
%       11          - 'one-third-OB'        YES         NO
%       12          - 'sharpness'           NO          NO
%       12          - 'loudness'            NO          NO
%       12          - 'specific-loudness'   YES         NO
%       15          - 'roughness'           NO          YES
%       15          - 'specific-roughness'  YES         YES
% 
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       PsySoundCL_Figures;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/08/2014
% Last update on: 01/10/2014 % Update this date manually
% Last use on   : 01/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
ha = [];
stats = [];

t   = res1.t; % assuming the same time vector

try 
    z   = res1.z;
end

try % not available for Roughness analyser yet...
    zspec = res1.zspec;
end

try
    f = res1.f;
end

if nargin < 5
    option = [];
end

stPlot = Ensure_field(stPlot,'label1','audio-1');
stPlot = Ensure_field(stPlot,'label2','audio-2');
stPlot = Ensure_field(stPlot,'label3','audio-3');
stPlot = Ensure_field(stPlot,'label1suffix','');
stPlot = Ensure_field(stPlot,'label2suffix','');
stPlot = Ensure_field(stPlot,'label3suffix','');

option = Ensure_field(option,'tanalysis',[min(res1.t) max(res1.t)]);

stPlot = Ensure_field(stPlot,'Title',[]);
option = Ensure_field(option,'color',{'b-','r--'});
option = Ensure_field(option,'LineWidth',[2 1]);

idx = find(t>=option.tanalysis(1) & t<=option.tanalysis(2));

if strcmp(param,'FFT')
    
    bPlot_vs_time = 0;
    Data1 = res1.Data1;
    Data2 = res2.Data1;
    freq_min = min(f);
    freq_max = max(f);
    
    idx = find(t>=option.tanalysis(1) & t<=option.tanalysis(2));

    figure;
    semilogx(f,mean(Data1(idx,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
    semilogx(f,mean(Data2(idx,:)),option.color{2},'LineWidth',option.LineWidth(2));
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)');
    strTitle = sprintf('Avg Power Spectrum - %s', stPlot.Title);
    if option.tanalysis(1) ~= 0
        strTitle = sprintf('%s\n, (ti,tf)=(%.3f,%.3f) [s]',strTitle,option.tanalysis(1),option.tanalysis(2))
    end
    title(strTitle);
    grid on;
    h(end+1) = gcf;
    ha(end+1) = gca;
    
elseif strcmp(param,'sharpness')
    
    bPlot_vs_time = 1;
    DataSharp1 = res1.DataSharp;
    DataSharp2 = res2.DataSharp;
        
    figure;
    plot(t,DataSharp1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t,DataSharp2,option.color{2},'LineWidth',option.LineWidth(2));

    xlabel('Time (Seconds)')
    ylabel('Sharpness (Acums)');
    title(sprintf('Sharpness - %s', stPlot.Title));
    grid on;
    h(end+1) = gcf;
    ha(end+1) = gca;

elseif strcmp(param,'loudness')
    % Loudness
    
    bPlot_vs_time = 1;
    DataLoud1 = res1.DataLoud;
    DataLoud2 = res2.DataLoud;
    DataLoud3 = res3.DataLoud;
    
    figure;
    plot(t,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
    plot(t,DataLoud3,option.color{3},'LineWidth',option.LineWidth(3));
    
    xlabel('Time (Seconds)')
    ylabel('Loudness (Sones)');
    title(sprintf('Loudness - %s', stPlot.Title));
    grid on;
    h(end+1) = gcf;
    ha(end+1) = gca;
 
elseif strcmp(param,'one-third-OB')
    
    bPlot_vs_time = 0;
    option = Ensure_field(option,'nAnalyser',10);
    
    switch option.nAnalyser
        case 10
            
            % nParam = 2;
            f = res1.f;
            freq_min = min(f);
            freq_max = max(f);
            DataSpecOneThirdAvg1 = res1.DataSpecOneThirdAvg;
            DataSpecOneThirdAvg2 = res2.DataSpecOneThirdAvg;
        
            figure;
            semilogx(f,DataSpecOneThirdAvg1,option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
            semilogx(f,DataSpecOneThirdAvg2,option.color{2},'LineWidth',option.LineWidth(2),'Marker','<'); grid on;
            ylabel('Magnitude (dB)')
            xlabel('Frequency (Hz)');
            title(sprintf('Average one-third octave band spectrum - %s', stPlot.Title));
            h(end+1) = gcf;
            ha(end+1) = gca;
            
        case 11
            
            f = res1.f;
            freq_min = min(f);
            freq_max = max(f);
            DataSpecOneThirdAvg1 = res1.DataSpecOneThirdAvg;
            DataSpecOneThirdAvg2 = res2.DataSpecOneThirdAvg;
        
            figure;
            
            if length(idx) == length(t)
                
                semilogx(f,DataSpecOneThirdAvg1,option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
                semilogx(f,DataSpecOneThirdAvg2,option.color{2},'LineWidth',option.LineWidth(2),'Marker','<'); grid on;
                title(sprintf('Average one-third octave band spectrum - %s', stPlot.Title));
                
            else
                
                title(sprintf('Average one-third octave band spectrum - %s (ti, tf) = (%.3f,%.3f) [s]', stPlot.Title,option.tanalysis(1),option.tanalysis(2)));   
                semilogx(f,dbmean( res1.DataSpecOneThird(idx,:) ),option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
                semilogx(f,dbmean( res2.DataSpecOneThird(idx,:) ),option.color{2},'LineWidth',option.LineWidth(2),'Marker','<'); grid on;
            
            end
            ylabel('Magnitude (dB)')
            xlabel('Frequency (Hz)');
                
            h(end+1) = gcf;
            ha(end+1) = gca;
            
    end
    
elseif strcmp(param,'specific-loudness')
    
    option = Ensure_field(option,'nAnalyser',12);
    bPlot_vs_time = 0;
    
    switch option.nAnalyser
        case 10
            
            stPlot.label1suffix = sprintf(', tot = %.2f [sone]',res1.stats.loud_tot);
            stPlot.label2suffix = sprintf(', tot = %.2f [sone]',res2.stats.loud_tot);

            % nParam = 3;
            zspec = res1.zspec;
            freq_min = min(zspec);
            freq_max = max(zspec);
            
            DataLoud1 = res1.DataLoud;
            DataLoud2 = res2.DataLoud;
        
            figure;
            plot(zspec,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
            plot(zspec,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
            xlabel('Critical band rate (Bark)')
            ylabel('Loudness (Sones/Bark)');
            
            title(sprintf('Specific Loudness (ISO532B) - %s', stPlot.Title));
            
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
            
        otherwise % analyser 12
            
            DataLoud1 = res1.DataAvSpecLoud;
            DataLoud2 = res2.DataAvSpecLoud;
            DataLoud3 = res3.DataAvSpecLoud;
            DataSpecLoud1 = res1.DataSpecLoud;
            DataSpecLoud2 = res2.DataSpecLoud;
            DataSpecLoud3 = res3.DataSpecLoud;
            
            figure;
            if length(idx) == length(t)
                
                plot(zspec,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
                plot(zspec,DataLoud3,option.color{3},'LineWidth',option.LineWidth(3));
                title(sprintf('Average Specific Loudness - %s', stPlot.Title));
                
            else
                % Excerpt
                plot(zspec,mean(DataSpecLoud1(idx,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,mean(DataSpecLoud2(idx,:)),option.color{2},'LineWidth',option.LineWidth(2));
                plot(zspec,mean(DataSpecLoud3(idx,:)),option.color{3},'LineWidth',option.LineWidth(3));
                title(sprintf('Average Specific Loudness - %s (ti, tf) = (%.3f,%.3f) [s]', stPlot.Title,option.tanalysis(1),option.tanalysis(2)));
                
            end
            
            xlabel('Critical band rate (Bark)');
            ylabel('Loudness (Sone/Bark)')
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
            
    end
    
elseif strcmp(param,'roughness')
    
    bPlot_vs_time = 1;
    DataRough1 = res1.DataRough;
    DataRough2 = res2.DataRough;

    figure;
    plot(t,DataRough1,option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
    plot(t,DataRough2,option.color{2},'LineWidth',option.LineWidth(2),'Marker','<');

    xlabel('Time (seconds)')
    ylabel('Roughness (aspers)')
    title(sprintf('Roughness - %s', stPlot.Title));
    grid on

    h(end+1) = gcf;
    ha(end+1) = gca;
        
elseif strcmp(param,'specific-roughness')
    
    bPlot_vs_time = 0;
    freq_min = min(z);
    freq_max = max(z);
    
    DataRough1 = res1.DataSpecRough;
    DataRough2 = res2.DataSpecRough;
    
    res1.stats.rough_segment = 0.5*sum(mean(DataRough1(idx,:)));
    res2.stats.rough_segment = 0.5*sum(mean(DataRough2(idx,:)));
    
    figure;
    plot(z, mean(DataRough1(idx,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z, mean(DataRough2(idx,:)),option.color{2},'LineWidth',option.LineWidth(2));
    xlabel('Critical band rate (Bark)')
    ylabel('Specific Roughness (Aspers/Bark)')
    %title(sprintf('Average Roughness - %s', stPlot.Title));
    
    if length(idx) == length(t)
        title(sprintf('Average Roughness - %s', stPlot.Title));
        stPlot.label1suffix = sprintf(', tot = %.2f [asper]',res1.stats.rough_tot);
        stPlot.label2suffix = sprintf(', tot = %.2f [asper]',res2.stats.rough_tot);
    else
        title(sprintf('Average Roughness - %s (ti, tf) = (%.3f,%.3f) [s]', stPlot.Title,option.tanalysis(1),option.tanalysis(2)));
        stPlot.label1suffix = sprintf(', tot = %.2f [asper]',res1.stats.rough_segment);
        stPlot.label2suffix = sprintf(', tot = %.2f [asper]',res2.stats.rough_segment);
    end
    
    grid on
    h(end+1) = gcf;
    ha(end+1) = gca;
    
    disp('sum(DataRough1): ');
    disp(['  : ' num2str(sum(mean(DataRough1))*0.1) ]);
    disp('sum(DataRough2): ');
    disp(['  : ' num2str(sum(mean(DataRough2))*0.1) ]);
    disp('sum(DataRough1-DataRough2): ');
    disp(['  : ' num2str(sum(mean(DataRough1)-mean(DataRough2))*0.1) ]);
    
end

stats.idx   = idx;
stats.t     = t;

legend( sprintf('%s %s',stPlot.label1,stPlot.label1suffix), ... 
        sprintf('%s %s',stPlot.label2,stPlot.label2suffix), ...
        sprintf('%s %s',stPlot.label3,stPlot.label3suffix));
    
try
    if bPlot_vs_time == 0
        xlim([freq_min freq_max]) % Bark
    else
        xlim(option.tanalysis)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
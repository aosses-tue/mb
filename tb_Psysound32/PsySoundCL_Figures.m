function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% 
% 1. Description:
%       To plot PsySound results for 2 data series
%       'param' can be:
% 
%       nAnalyser                           excerpt     stats
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

if nargin < 4
    option = [];
end

option = Ensure_field(option,'label1','audio-1');
option = Ensure_field(option,'label2','audio-2');
option = Ensure_field(option,'label1suffix','');
option = Ensure_field(option,'label2suffix','');

option = Ensure_field(option,'tanalysis',[min(res1.t) max(res1.t)]);

option = Ensure_field(option,'title',[]);
option = Ensure_field(option,'color',{'b-','r--'});
option = Ensure_field(option,'LineWidth',[2 1]);

idx = find(t>=option.tanalysis(1) & t<=option.tanalysis(2));

if strcmp(param,'sharpness')
    
    bPlot_vs_time = 1;
    DataSharp1 = res1.DataSharp;
    DataSharp2 = res2.DataSharp;
        
    figure;
    plot(t,DataSharp1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t,DataSharp2,option.color{2},'LineWidth',option.LineWidth(2));

    xlabel('Time (Seconds)')
    ylabel('Sharpness (Acums)');
    title(sprintf('Sharpness - %s', option.title));
    grid on;
    h(end+1) = gcf;
    ha(end+1) = gca;

elseif strcmp(param,'loudness')
    % Loudness
    
    bPlot_vs_time = 1;
    DataLoud1 = res1.DataLoud;
    DataLoud2 = res2.DataLoud;

    figure;
    plot(t,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));

    xlabel('Time (Seconds)')
    ylabel('Loudness (Sones)');
    title(sprintf('Loudness - %s', option.title));
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
            title(sprintf('Average one-third octave band spectrum - %s', option.title));
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
                title(sprintf('Average one-third octave band spectrum - %s', option.title));
                
            else
                
                title(sprintf('Average one-third octave band spectrum - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));   
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
            
            option.label1suffix = sprintf(', tot = %.2f [sone]',res1.stats.loud_tot);
            option.label2suffix = sprintf(', tot = %.2f [sone]',res2.stats.loud_tot);

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
            
            title(sprintf('Specific Loudness (ISO532B) - %s', option.title));
            
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
            
        otherwise % analyser 12
            
            zspec = res1.zspec;
            freq_min = min(zspec);
            freq_max = max(zspec);
            
            DataLoud1 = res1.DataAvSpecLoud;
            DataLoud2 = res2.DataAvSpecLoud;
            DataSpecLoud1 = res1.DataSpecLoud;
            DataSpecLoud2 = res2.DataSpecLoud;
            
            figure;
            if length(idx) == length(t)
                
                plot(zspec,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
                title(sprintf('Average Specific Loudness - %s', option.title));
                data2show1 = sum(DataLoud1)*0.1;
                data2show2 = sum(DataLoud2)*0.1;
                
            else
                % Excerpt
                plot(zspec,mean(DataSpecLoud1(idx,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,mean(DataSpecLoud2(idx,:)),option.color{2},'LineWidth',option.LineWidth(2));
                title(sprintf('Average Specific Loudness - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
                data2show1 = sum(DataSpecLoud1(idx,:))*0.1;
                data2show2 = sum(DataSpecLoud2(idx,:))*0.1;
                
            end
            
            xlabel('Critical band rate (Bark)');
            ylabel('Loudness (Sone/Bark)')
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
            
            disp('sum(DataLoud1): ');
            disp(['  : ' num2str(data2show1) ]);
            disp('sum(DataLoud2): ');
            disp(['  : ' num2str(data2show2) ]);
            disp('sum(DataLoud1-DataLoud2): ');
            disp(['  : ' num2str(data2show1 - data2show2) ]);
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
    title(sprintf('Roughness - %s', option.title));
    grid on

    res1.stats.rough_segment = mean(DataRough1(idx));
    res2.stats.rough_segment = mean(DataRough2(idx));
    
    h(end+1) = gcf;
    ha(end+1) = gca;
        
elseif strcmp(param,'specific-roughness')
    
    bPlot_vs_time = 0;
    freq_min = min(z);
    freq_max = max(z);
    
    DataRough1 = res1.DataSpecRough;
    DataRough2 = res2.DataSpecRough;
    
    res1.stats.rough_segment = mean( 0.25*sum( DataRough1(idx,:)' ) ); % 0.5*sum( mean(DataRough1(idx,:)') );
    res2.stats.rough_segment = mean( 0.25*sum( DataRough2(idx,:)' ) );
    
    figure;
    plot(z, mean(DataRough1(idx,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z, mean(DataRough2(idx,:)),option.color{2},'LineWidth',option.LineWidth(2));
    xlabel('Critical band rate (Bark)')
    ylabel('Specific Roughness (Aspers/Bark)')
    %title(sprintf('Average Roughness - %s', option.title));
    
    if length(idx) == length(t)
        title(sprintf('Average Roughness - %s', option.title));
        option.label1suffix = sprintf(', tot = %.2f [asper]',res1.stats.rough_tot);
        option.label2suffix = sprintf(', tot = %.2f [asper]',res2.stats.rough_tot);
    else
        title(sprintf('Average Roughness - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
        option.label1suffix = sprintf(', tot = %.2f [asper]',res1.stats.rough_segment);
        option.label2suffix = sprintf(', tot = %.2f [asper]',res2.stats.rough_segment);
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

legend( sprintf('%s %s',option.label1,option.label1suffix), ... 
        sprintf('%s %s',option.label2,option.label2suffix) );
    
try
    if bPlot_vs_time == 0
        xlim([freq_min freq_max]) % Bark
    else
        xlim(option.tanalysis)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
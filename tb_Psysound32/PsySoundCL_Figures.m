function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% 
% 1. Description:
%       To plot PsySound results for 2 data series
%       'param' can be:
% 
%       nAnalyser                           excerpt     stats
%       1           - FFT                   PsySoundCL_Figures not working, probably one of the revisions was overwritten (see r20141007_Perception_day_TUe.m)
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
% Last use on   : 21/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
ha = [];
stats = [];

try
    t1   = res1.t; 
    t2   = res2.t; 
catch
    warning('no time variable was found...')
end

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
option = ef(option,'bGenerateExcerpt',0);

minValue = max( min(res1.t), min(res2.t) ); % normally = 0
maxValue = min( max(res1.t), max(res2.t) );

if length(res1.t) ~= length(res2.t)
    warning('Average values based on truncated time series...');
end

option = Ensure_field(option,'tanalysis',[minValue maxValue]);

option = Ensure_field(option,'title',[]);
option = Ensure_field(option,'color',{'b-','r--'});
option = Ensure_field(option,'LineWidth',[2 1]);

idx1 = find(t1>=option.tanalysis(1) & t1<=option.tanalysis(2));
idx2 = find(t2>=option.tanalysis(1) & t2<=option.tanalysis(2));

if option.bGenerateExcerpt == 1
    idx1 = find(t1>=0 & t1<=option.tanalysis(2)-option.tanalysis(1));
    idx2 = find(t2>=0 & t2<=option.tanalysis(2)-option.tanalysis(1));
    timeoffset = option.tanalysis(1);
else
    timeoffset = 0;
end

if strcmp(param,'sharpness')
    
    bPlot_vs_time = 1;
    DataSharp1 = res1.DataSharp;
    DataSharp2 = res2.DataSharp;
        
    figure;
    plot(t1,DataSharp1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t2,DataSharp2,option.color{2},'LineWidth',option.LineWidth(2));

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
    plot(t1,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t2,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));

    xlabel('Time (Seconds)')
    ylabel('Loudness (Sones)');
    title(sprintf('Loudness - %s', option.title));
    grid on;
    h(end+1) = gcf;
    ha(end+1) = gca;
 
elseif strcmp(param,'loudness-fluctuation')
    
    bPlot_vs_time = 0;
    Data1max = res1.Data1;
    Data2max = res2.Data1;
    Data1min = res1.Data2;
    Data2min = res2.Data2;
    
    figure;
    subplot(2,2,1)
    plot(z,Data1max,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z,Data2max,option.color{2},'LineWidth',option.LineWidth(2));
    grid on
    legend(option.label1,option.label2)
    
    xlabel('Critical band rate (Bark)')
    ylabel('Critical-band level L_G [dB]')
    ha = gca;
    title('L_G_{max}')
    
    subplot(2,2,2)
    plot(z,Data1min,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z,Data2min,option.color{2},'LineWidth',option.LineWidth(2));

    grid on
    
    xlabel('Critical band rate (Bark)')
    ylabel('Critical-band level L_G [dB]')
    % title(sprintf('Loudness - %s', option.title));
    % ylim(stPlot.YLim_fig1)
    % title(stPlot.Title2)
    
    ha(end+1) = gca;
    title('L_G_{min}')
    
    % title('L_G_{max} based on N_{95}')
    
    subplot(2,2,3)
    plot(   z,Data1max-Data2max,'bo-')
    ylabel('\Delta L_G [dB]')
    xlabel('Critical-band rate [Bark]')
    grid on
    % ylim(stPlot.YLim_fig2)
    % title(stPlot.Title4)
    ha(end+1) = gca;
    legend('diff max. levels')
    
    subplot(2,2,4)
    plot(   z,Data1min-Data2min,'bo-')
    ylabel('\Delta L_G [dB]')
    xlabel('Critical-band rate [Bark]')
    grid on
    
    % ylim(stPlot.YLim_fig2)
    % title(stPlot.Title4)
    h(end+1) = gcf;
    ha(end+1) = gca;
    legend('diff min. levels')
    
    linkaxes(ha,'x');
    % xlim([0 24])
    
elseif strcmp(param,'one-third-OB')
    
    bPlot_vs_time = 0;
    option = Ensure_field(option,'nAnalyser',10);
    
    switch option.nAnalyser
            
        case 10
            
            % nParam = 2;
            f = res1.f;
            option = Ensure_field(option,'frange',minmax(f));
            freq_min = min(f);
            freq_max = max(f);
            DataSpecOneThirdAvg1 = res1.DataSpecOneThirdAvg;
            DataSpecOneThirdAvg2 = res2.DataSpecOneThirdAvg;
        
            figure;
            semilogx(f,DataSpecOneThirdAvg1,option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
            semilogx(f,DataSpecOneThirdAvg2,option.color{2},'LineWidth',option.LineWidth(2),'Marker','<'); grid on;
            
            xlim(option.frange);
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
            
            if length(idx1) == length(t1)
                
                semilogx(f,DataSpecOneThirdAvg1,option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
                semilogx(f,DataSpecOneThirdAvg2,option.color{2},'LineWidth',option.LineWidth(2),'Marker','<'); grid on;
                title(sprintf('Average one-third octave band spectrum - %s', option.title));
                
            else
                
                title(sprintf('Average one-third octave band spectrum - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));   
                semilogx(f,dbmean( res1.DataSpecOneThird(idx1,:) ),option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
                semilogx(f,dbmean( res2.DataSpecOneThird(idx1,:) ),option.color{2},'LineWidth',option.LineWidth(2),'Marker','<'); grid on;
            
            end
            ylabel('Magnitude (dB)')
            xlabel('Frequency (Hz)');
                
            h(end+1) = gcf;
            ha(end+1) = gca;
            
    end
    
elseif strcmp(param,'average-power-spectrum')
    
    bPlot_vs_time = 0;
    % option = Ensure_field(option,'nAnalyser',10);
    
    switch option.nAnalyser
        case 1
            f = res1.f;
            option = ef(option,'frange',minmax(f));
            
            figure;
            semilogx(f,dbmean( res1.Data1(:,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
            semilogx(f,dbmean( res2.Data1(:,:)),option.color{2},'LineWidth',option.LineWidth(2)); grid on;
            ylabel('Magnitude (dB)'); % this can be automised
            xlabel('Frequency (Hz)'); % this can be automised
            xlim(option.frange)
            
            title(sprintf('%s - %s (ti, tf) = (%.3f,%.3f) [s]',res1.name{2},option.title,option.tanalysis(1),option.tanalysis(2)));
            h(end+1) = gcf;
            ha(end+1) = gca;
    end 
    
elseif strcmp(param,'spectrogram')
    
    switch option.nAnalyser
        case 1
            f = res1.f;
            option = ef(option,'frange',minmax(f));
            t = transpose( res1.t );
            
            figure;
            subplot(1,2,1)
            imagesc( t(idx1), f, transpose(res1.Data1(idx1,:)) );
            colormap('Gray')
            set(gca,'YDir','normal')
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
            xlabel('Time (s)'); 
            ylabel('Frequency (Hz)'); 
            h(end+1) = gcf;
            ha(end+1) = gca;
            title( sprintf('%s %s',option.label1,option.label1suffix) );
            ylim(option.frange);
            
            subplot(1,2,2)
            imagesc( t(idx2), f, transpose(res2.Data1(idx2,:)) );
            colormap('Gray')
            set(gca,'YDir','normal')
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
            xlabel('Time (s)'); 
            ylabel('Frequency (Hz)'); 
            ylim(option.frange);
            % h(end+1) = gcf;
            % ha(end+1) = gca;
            title( sprintf('%s %s',option.label2,option.label2suffix) );
            
            linkaxes([ha gca],'xy');
    end
    
elseif strcmp(param,'specific-loudness')| strcmp(param,'average-specific-loudness')
    
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
            if length(idx1) == length(t1)
                
                plot(zspec,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
                title(sprintf('Average Specific Loudness - %s', option.title));
                data2show1 = sum(DataLoud1)*0.1;
                data2show2 = sum(DataLoud2)*0.1;
                
            else
                % Excerpt
                plot(zspec,mean(DataSpecLoud1(idx1,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,mean(DataSpecLoud2(idx1,:)),option.color{2},'LineWidth',option.LineWidth(2));
                title(sprintf('Average Specific Loudness - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
                data2show1 = sum(DataSpecLoud1(idx1,:))*0.1;
                data2show2 = sum(DataSpecLoud2(idx1,:))*0.1;
                
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
    Data1 = res1.Data1; % res1.DataRough;
    Data2 = res2.Data1; % res2.DataRough;
    
    figure;
    plot(t1+timeoffset,Data1,option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
    plot(t2+timeoffset,Data2,option.color{2},'LineWidth',option.LineWidth(2),'Marker','<');

    xlabel('Time (seconds)')
    ylabel('Roughness (aspers)')
    title(sprintf('Roughness - %s', option.title));
    grid on

    res1.stats.rough_segment = mean(Data1(idx1));
    res2.stats.rough_segment = mean(Data2(idx2));
    
    h(end+1) = gcf;
    ha(end+1) = gca;
        
elseif strcmp(param,'specific-roughness')| strcmp(param,'average-specific-roughness')
    
    bPlot_vs_time = 0;
    freq_min = min(z);
    freq_max = max(z);
    
    Data1 = res1.Data2; % res1.DataSpecRough;
    Data2 = res2.Data2; % res2.DataSpecRough;
    
    res1.stats.rough_segment = mean( 0.25*sum( Data1(idx1,:)' ) ); % 0.5*sum( mean(DataRough1(idx,:)') );
    res2.stats.rough_segment = mean( 0.25*sum( Data2(idx1,:)' ) );
    
    figure;
    plot(z, mean(Data1(idx1,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z, mean(Data2(idx2,:)),option.color{2},'LineWidth',option.LineWidth(2));
    xlabel('Critical band rate (Bark)')
    ylabel('Specific Roughness (Aspers/Bark)')
    %title(sprintf('Average Roughness - %s', option.title));
    
    if length(idx1) == length(t1)
        try
            title(sprintf('Average Roughness - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
        catch
            title(sprintf('Average Roughness - %s', option.title));
        end
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
    
end

stats.idx   = idx1;
stats.t     = t1;

if ~strcmp(param,'spectrogram')
    legend( sprintf('%s %s',option.label1,option.label1suffix), ... 
            sprintf('%s %s',option.label2,option.label2suffix) );
end
    
try
    if bPlot_vs_time == 0
        xlim([freq_min freq_max]) % Bark
    else
        % xlim(option.tanalysis)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
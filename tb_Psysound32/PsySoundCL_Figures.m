function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% 
% 1. Description:
%       To plot PsySound results for 2 data series
%       'param' can be:
% 
%       nAnalyser                           excerpt     stats
% (XX/XX/2015)      1       - FFT                   PsySoundCL_Figures not working, probably one of the revisions was overwritten (see r20141007_Perception_day_TUe.m)
% (XX/XX/2015)     10       - 'one-third-OB'        NO          NO
% (XX/XX/2015)     10       - 'specific-loudness'   NO          NO
% (XX/XX/2015)        11    - 'one-third-OB'        YES         NO
% (25/02/2015) L96  / 12    - 'sharpness'           NO          NO
% (25/02/2015) L117 / 12    - 'loudness'            NO          NO
% (25/02/2015) L417 / 12    - 'loudness-percentiles'
% (01/03/2015) L345 / 12    - 'specific-loudness'   YES         NO
% (05/02/2015) L475 / 15    - 'roughness'           NO          YES
%       15                  - 'specific-roughness'  YES         YES
% (29/01/2015) L249 /       - 'average-power-spectrum' 
% (25/02/2015) L310 /       - 'spectrogram'
% 
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       PsySoundCL_Figures;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/08/2014
% Last update on: 16/02/2015 % Update this date manually
% Last use on   : 16/02/2015 % Update this date manually
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

option = Ensure_field(option,'zrange',[0 24]);

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
option = ef(option,'bUsePsySound'    ,1);
option = ef(option,'bLoudnessContrained',0);
option = ef(option,'zlim4assessment',[0 24]);

bLoudnessContrained = option.bLoudnessContrained; % out of the struct to be visible by MATLAB
zlim4assessment = option.zlim4assessment;

minValue = max( min(res1.t), min(res2.t) ); % normally = 0
maxValue = min( max(res1.t), max(res2.t) );

if length(res1.t) ~= length(res2.t)
    warning('Average values based on truncated time series...');
end

option = Ensure_field(option,'tanalysis',[minValue maxValue]);

option = ef(option,'trange',option.tanalysis);
option = Ensure_field(option,'title',[]);
option = Ensure_field(option,'color'          ,{'b-','r--'});
option = Ensure_field(option,'color_no_format',{'b' ,'r'}  );
option = Ensure_field(option,'LineWidth',[2 1]);

idx1 = find(t1>=option.tanalysis(1) & t1<=option.tanalysis(2));
idx2 = find(t2>=option.tanalysis(1) & t2<=option.tanalysis(2));

if option.bGenerateExcerpt == 1 & option.nAnalyser == 15 & option.bUsePsySound == 0 % only for roughness offline
    idx1 = find(t1>=0 & t1<=option.tanalysis(2)-option.tanalysis(1));
    idx2 = find(t2>=0 & t2<=option.tanalysis(2)-option.tanalysis(1));
     timeoffset = option.tanalysis(1);
else
    timeoffset = 0;
end

if strcmp(param,'sharpness')
    
    bPlot_vs_time = 1;
    DataSharp1 = res1.Data6;
    DataSharp2 = res2.Data6;
        
    figure;
    plot(t1,DataSharp1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t2,DataSharp2,option.color{2},'LineWidth',option.LineWidth(2));

    xlabel('Time (Seconds)')
    ylabel('Sharpness (Acums)');
    if bLoudnessContrained == 0
        title(sprintf('Sharpness - %s', option.title));
    else
        title(sprintf('Sharpness, z = (%.1f,%.1f) [Bark] - %s', zlim4assessment(1),zlim4assessment(2),option.title));
    end
    grid on;
    h(end+1) = gcf;
    ha(end+1) = gca;

elseif strcmp(param,'loudness')
    % Loudness
    
    bPlot_vs_time = 1;
    DataLoud1 = res1.Data1;
    DataLoud2 = res2.Data1;

    figure;
    plot(t1,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t2,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));

    xlabel('Time (Seconds)')
    ylabel('Loudness (Sones)');
    if bLoudnessContrained == 0
        title(sprintf('Loudness - %s', option.title));
    else
        title(sprintf('Loudness, z = (%.1f,%.1f) [Bark] - %s', zlim4assessment(1),zlim4assessment(2),option.title));
    end
    grid on;
    
    h(end+1) = gcf;
    ha(end+1) = gca;
 
elseif strcmp(param,'loudness-fluctuation')| strcmp(param,'loudness-fluctuation-max')
    
    bPlot_vs_time = 0;
    Data1max = res1.Data1;
    Data2max = res2.Data1;
    
    figure;
    subplot(2,1,1)
    plot(z,Data1max,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z,Data2max,option.color{2},'LineWidth',option.LineWidth(2));
    grid on
    legend(option.label1,option.label2)
    
    xlabel('Critical band rate (Bark)')
    ylabel('Critical-band level L_G [dB]')
    ha = gca;
    title('L_G_{max}')
    
    subplot(2,1,2)
    plot(   z,Data1max-Data2max,'bo-')
    ylabel('\Delta L_G [dB]')
    xlabel('Critical-band rate [Bark]')
    grid on
    legend('diff max. levels')
    
    ha(end+1) = gca;
    h(end+1) = gcf;
    
    linkaxes(ha,'x');
    xlim(option.zrange);
    
elseif strcmp(param,'loudness-fluctuation-min')
    
    bPlot_vs_time = 0;
    Data1min = res1.Data2;
    Data2min = res2.Data2;
    
    figure
    subplot(2,1,1)
    plot(z,Data1min,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z,Data2min,option.color{2},'LineWidth',option.LineWidth(2));

    grid on
    legend(option.label1,option.label2)
    
    xlabel('Critical band rate (Bark)')
    ylabel('Critical-band level L_G [dB]')
    ha(end+1) = gca;
    title('L_G_{min}')
    
    subplot(2,1,2)
    plot(   z,Data1min-Data2min,'bo-')
    ylabel('\Delta L_G [dB]')
    xlabel('Critical-band rate [Bark]')
    grid on
    legend('diff min. levels')
    
    ha(end+1) = gca;
    h(end+1) = gcf;
        
    linkaxes(ha,'x');
    xlim(option.zrange);
    
elseif strcmp(param,'one-third-OB') | strcmp(param,'one-third-octave-band-spectrum')
    
    bPlot_vs_time = 0;
    option = Ensure_field(option,'nAnalyser',10);
    
    switch option.nAnalyser
            
        case 10
            
            % nParam = 2;
            f = res1.f;
            option = Ensure_field(option,'frange',minmax(f));
            freq_min = option.frange(1);
            freq_max = option.frange(2);
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
                
            ylims = get(gca,'YLim');
            ylims(1) = 0; % Nothing below 0 dB
            set(gca,'YLim',ylims);
            
            xticks = get(gca,'XTick');
            set(gca,'XTickLabel',xticks);
            
            h(end+1) = gcf;
            ha(end+1) = gca;
            
    end
    
elseif strcmp(param,'average-power-spectrum')
    
    bPlot_vs_time = 0;
    % option = Ensure_field(option,'nAnalyser',10);
    
    switch option.nAnalyser
        case 1
            option = Ensure_field(option,'bLogScale',1);
            f = res1.f;
            option = ef(option,'frange',minmax(f));
            
            figure;
            if option.bLogScale == 1
                semilogx(f,dbmean( res1.Data1(:,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
                semilogx(f,dbmean( res2.Data1(:,:)),option.color{2},'LineWidth',option.LineWidth(2)); grid on;
            else
                plot(f,dbmean( res1.Data1(:,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(f,dbmean( res2.Data1(:,:)),option.color{2},'LineWidth',option.LineWidth(2)); grid on;
            end
            ylabel('Magnitude (dB)'); % this can be automised
            xlabel('Frequency (Hz)'); % this can be automised
            xlim(option.frange)
            
            title(sprintf('%s - %s (ti, tf) = (%.3f,%.3f) [s]',res1.name{2},option.title,option.tanalysis(1),option.tanalysis(2)));
            
            ylims = get(gca,'YLim');
            ylims(1) = 0; % Nothing below 0 dB
            set(gca,'YLim',ylims);
            
            xticks = get(gca,'XTick');
            % set(gca,'XTickLabel',xticks);
            
            hhAxes = handle(gca);  % hAxes is the Matlab handle of our axes
            hProp = findprop(hhAxes,'XTick');  % a schema.prop object
            hListener = handle.listener(hhAxes, hProp, 'PropertyPostSet', @myCallbackFunction);
            setappdata(gca, 'XTickListener', hListener);
            
            h(end+1) = gcf;
            ha(end+1) = gca;
    end 
    
elseif strcmp(param,'spectrogram')
    bPlot_vs_time = 1;
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
            title( sprintf('%s %s',option.label2,option.label2suffix) );
            
            linkaxes([ha gca],'xy');
    end
    
elseif strcmp(param,'specific-loudness')| strcmp(param,'average-specific-loudness')
    
    option = Ensure_field(option,'nAnalyser',12);
        
    switch option.nAnalyser
        case 10
            
            option.label1suffix = sprintf(', tot = %.2f [sone]',res1.stats.loud_tot);
            option.label2suffix = sprintf(', tot = %.2f [sone]',res2.stats.loud_tot);

            % nParam = 3;
            zspec = res1.zspec;
            freq_min = min(zspec);
            freq_max = max(zspec);
            
            DataLoud1 = res1.Data1;
            DataLoud2 = res2.Data1;
        
            figure;
            plot(zspec,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
            plot(zspec,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
            xlabel('Critical band rate (Bark)')
            ylabel('Loudness (Sones/Bark)');
            xlim(option.zrange);
            
            title(sprintf('Specific Loudness (ISO532B) - %s', option.title));
            
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;
            
        otherwise % analyser 12
            
            zspec = res1.zspec;
            freq_min = min(zspec);
            freq_max = max(zspec);
            
            DataLoud1 = res1.Data5;
            DataLoud2 = res2.Data5;
            DataSpecLoud1 = res1.Data3; % Spec-loud
            DataSpecLoud2 = res2.Data3; % Spec-loud
            
            figure;
            if length(idx1) == length(t1)
                
                plot(zspec,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
                title(sprintf('Average Specific Loudness - %s', option.title));
                data2show1 = sum(DataLoud1)*0.1;
                data2show2 = sum(DataLoud2)*0.1;
                xlim(option.zrange);
                
            else
                % Excerpt
                plot(zspec,mean(DataSpecLoud1(idx1,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
                plot(zspec,mean(DataSpecLoud2(idx1,:)),option.color{2},'LineWidth',option.LineWidth(2));
                title(sprintf('Average Specific Loudness - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
                data2show1 = sum( mean(DataSpecLoud1(idx1,:)) )*0.1;
                data2show2 = sum( mean(DataSpecLoud2(idx1,:)) )*0.1;
                xlim(option.zrange);
                
            end
            
            option.label1suffix = sprintf(', avg=%.2f [sones]',data2show1);
            option.label2suffix = sprintf(', avg=%.2f [sones]',data2show2);
                
            xlabel('Critical band rate (Bark)');
            ylabel('Loudness (Sone/Bark)')
            grid on;
            h(end+1) = gcf;
            ha(end+1) = gca;

    end
elseif strcmp(param,'loudness-percentiles')
    
    idx = find(res1.t(idx1) >= option.trange(1) & res1.t(idx1) <= option.trange(2));
    
    N1 = 95;
    N2 = 50;
    N3 = 5;
    nmax1    = ch_prctile(res1.Data3(idx,:),N1);
    nmean1   = ch_prctile(res1.Data3(idx,:),N2);
    nmin1    = ch_prctile(res1.Data3(idx,:),N3);

    nmax2    = ch_prctile(res2.Data3(idx,:),N1);
    nmean2   = ch_prctile(res2.Data3(idx,:),N2);
    nmin2    = ch_prctile(res2.Data3(idx,:),N3);
             
    figure
    subplot(2,1,1)
    plot(zspec,nmax1,option.color_no_format{1},'LineStyle','--'), hold on %,'Marker','>')
    plot(zspec,nmean1,option.color_no_format{1},'LineWidth',2)
    plot(zspec,nmin1,option.color_no_format{1},'LineStyle','-.')%,'Marker','+')
    grid on
    
    xlabel('Critical-band rate (Bark)')
    ylabel('Specific loudness (Sones/Bark)')
    
    if bLoudnessContrained == 0
        title( sprintf('Percentiles %.0f, %.0f, %0.f, t =(%.1f,%.1f) [s] - %s %s',N3,N2,N1,option.trange(1),option.trange(2),option.label1,option.label1suffix));
    else
        title( sprintf('Percentiles %.0f, %.0f, %0.f; z = (%.1f,%.1f) [Bark], t =(%.1f,%.1f) [s] - %s %s',N3,N2,N1,zlim4assessment(1),zlim4assessment(2),option.trange(1),option.trange(2),option.label1,option.label1suffix));
    end
    
    legend( sprintf('N%.0f=%.1f [sone]',N1,0.1*sum(nmax1) ), ...
            sprintf('N%.0f=%.1f [sone]',N2,0.1*sum(nmean1)), ...
            sprintf('N%.0f=%.1f [sone]',N3,0.1*sum(nmin1) ) );
    
    grid on
    xlim(option.zrange);
    ha(end+1) = gca;
    
    subplot(2,1,2)
    plot(zspec,nmax2,option.color_no_format{2},'LineStyle','--'), hold on %,'Marker','>')
    plot(zspec,nmean2,option.color_no_format{2},'LineWidth',2)
    plot(zspec,nmin2,option.color_no_format{2},'LineStyle','-.')%,'Marker','+')
    grid on            
    
    xlabel('Critical-band rate (Bark)')
    ylabel('Specific loudness (Sones/Bark)')
    title( sprintf('Percentiles %.0f, %.0f, %0.f - %s %s',N3,N2,N1,option.label2,option.label2suffix));
    xlim(option.zrange);
    
    legend( sprintf('N%.0f=%.1f [sone]',N1,0.1*sum(nmax2)), ...
            sprintf('N%.0f=%.1f [sone]',N2,0.1*sum(nmean2)), ...
            sprintf('N%.0f=%.1f [sone]',N3,0.1*sum(nmin2)) );
    
    h(end+1) = gcf;
    ha(end+1) = gca;
    linkaxes(ha,'xy');
    
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
    xlim(option.zrange);
    
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

elseif strcmp(param,'fluctuation-strength')
    
    bPlot_vs_time = 1;
    Data1 = res1.Data1; 
    Data2 = res2.Data1; 
    
    figure;
    plot(t1+timeoffset,Data1,option.color{1},'LineWidth',option.LineWidth(1),'Marker','o'); hold on
    plot(t2+timeoffset,Data2,option.color{2},'LineWidth',option.LineWidth(2),'Marker','<');

    xlabel('Time (seconds)')
    ylabel('Fluctuation strength (vacils)')
    title(sprintf('FS - %s', option.title));
    grid on
    
    h(end+1) = gcf;
    ha(end+1) = gca;
        
elseif strcmp(param,'specific-fluctuation-strength')| strcmp(param,'average-specific-fluctuation-strength')
    
    bPlot_vs_time = 0;
    freq_min = min(z);
    freq_max = max(z);
    
    Data1 = res1.Data2; 
    Data2 = res2.Data2; 
    
    figure;
    if size(Data1,1) ~= 1
        plot(z, mean(Data1(idx1,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
        plot(z, mean(Data2(idx2,:)),option.color{2},'LineWidth',option.LineWidth(2));
    else % average of one value is the same value
        plot(z, Data1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
        plot(z, Data2,option.color{2},'LineWidth',option.LineWidth(2));
    end
    xlabel('Critical band rate (Bark)')
    ylabel('Specific Fluctuation strength (Vacils/Bark)')
    %title(sprintf('Average Roughness - %s', option.title));
    xlim(option.zrange);
    
    if length(idx1) == length(t1)
        try
            title(sprintf('Average FS - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
        catch
            title(sprintf('Average FS - %s', option.title));
        end
    else
        title(sprintf('Average FS - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));

    end
    
    grid on
    h(end+1) = gcf;
    ha(end+1) = gca;

end

stats.idx   = idx1;
stats.t     = t1;

if ~strcmp(param,'spectrogram') & ~strcmp(param,'loudness-fluctuation-max') & ... 
   ~strcmp(param,'loudness-fluctuation-min') & ~strcmp(param,'loudness-percentiles')
    legend( sprintf('%s %s',option.label1,option.label1suffix), ... 
            sprintf('%s %s',option.label2,option.label2suffix) );
end
    
try
    if bPlot_vs_time == 0
        xlim([freq_min freq_max]); % Bark
    else
        xlim(option.trange);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function myCallbackFunction(hProp,eventData)    %#ok - hProp is unused
   hAxes = eventData.AffectedObject;
   tickValues = get(hAxes,'XTick');
   newLabels = arrayfun(@(value)(sprintf('%.0f',value)), tickValues, 'UniformOutput',false);
   set(hAxes, 'XTickLabel', newLabels);
end  % myCallbackFunction
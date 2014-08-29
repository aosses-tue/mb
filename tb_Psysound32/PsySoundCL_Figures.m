function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% function [h ha stats] = PsySoundCL_Figures(param,res1, res2, option)
% 
% 1. Description:
%       To plot PsySound results for 2 data series
%       'param' can be:
%           - 'sharpness'
%           - 'loudness'
%           - 'roughness'
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       PsySoundCL_Figures;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/08/2014
% Last update on: 20/08/2014 % Update this date manually
% Last use on   : 29/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
ha = [];

t   = res1.t; % assuming the same time vector
z   = res1.z;

try % not available for Roughness analyser yet...
    zspec = res1.zspec;
end

if nargin < 4
    option = [];
end

option = Ensure_field(option,'tanalysis',[min(res1.t) max(res1.t)]);
option = Ensure_field(option,'title',[]);
option.color{1} = 'b-';
option.color{2} = 'r-';
option = Ensure_field(option,'LineWidth',[1 2]);

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
    stats = [];
        
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
    stats = [];

elseif strcmp(param,'specific-loudness')
    
    bPlot_vs_time = 0;
    DataLoud1 = res1.DataAvSpecLoud;
    DataLoud2 = res2.DataAvSpecLoud;
    
    figure;
    idx = find(t>=option.tanalysis(1) & t<=option.tanalysis(2));
    plot(zspec,DataLoud1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(zspec,DataLoud2,option.color{2},'LineWidth',option.LineWidth(2));
    xlabel('Critical band rate (Bark)');
    ylabel('Loudness (Sone/Bark)')
    
    if length(idx) == length(t)
        title(sprintf('Average Specific Loudness - %s', option.title));
    else
        title(sprintf('Average Specific Loudness - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
    end
    grid on;
    h(end+1) = gcf;
    ha(end+1) = gca;
    stats = [];
    
elseif strcmp(param,'roughness')
    
    bPlot_vs_time = 1;
    DataRough1 = res1.DataRough;
    DataRough2 = res2.DataRough;

    figure;
    plot(t,DataRough1,option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(t,DataRough2,option.color{2},'LineWidth',option.LineWidth(2));

    xlabel('Time (seconds)')
    ylabel('Roughness (aspers)')
    title(sprintf('Roughness - %s', option.title));
    grid on

    h(end+1) = gcf;
    ha(end+1) = gca;
    stats = [];
        
elseif strcmp(param,'specific-roughness')
    
    bPlot_vs_time = 0;
    DataRough1 = res1.DataSpecRough;
    DataRough2 = res2.DataSpecRough;
    
    figure;
    idx = find(t>=option.tanalysis(1) & t<=option.tanalysis(2));
    plot(z, mean(DataRough1(idx,:)),option.color{1},'LineWidth',option.LineWidth(1)); hold on
    plot(z, mean(DataRough2(idx,:)),option.color{2},'LineWidth',option.LineWidth(2));
    xlabel('Critical band rate (Bark)')
    ylabel('Specific Roughness (Aspers/Bark)')
    %title(sprintf('Average Roughness - %s', option.title));
    if length(idx) == length(t)
        title(sprintf('Average Roughness - %s', option.title));
    else
        title(sprintf('Average Roughness - %s (ti, tf) = (%.3f,%.3f) [s]', option.title,option.tanalysis(1),option.tanalysis(2)));
    end
    grid on
    h(end+1) = gcf;
    ha(end+1) = gca;
    stats = [];
    
end

try
    if bPlot_vs_time == 0
        xlim([0 24]) % Bark
    else
        xlim(option.tanalysis)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
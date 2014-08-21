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
% Last use on   : 20/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
ha = [];

t   = res1.t; % assuming the same time vector

if nargin < 4
    option = [];
end

option = Ensure_field(option,'tanalysis',[min(res1.t) max(res1.t)]);
option = Ensure_field(option,'title',[]);
option.color{1} = 'b-';
option.color{2} = 'r-';
option = Ensure_field(option,'LineWidth',[1 2]);

if strcmp(param,'sharpness')
        
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
        
elseif strcmp(param,'roughness')
    
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
end

try
    xlim(option.tanalysis)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
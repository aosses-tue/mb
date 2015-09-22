function r20150528_get_paper_plots
% function r20150528_get_paper_plots
%
% 1. Description:
%       Produces figures as added into the Euro-noise conference paper
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/02/2015
% Last update on: 28/05/2015 % Update this date manually
% Last use on   : 28/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

if isunix
    outputdir = '~/Documenten/Documenten-TUe/01-Text/05-Doc-TUe/lx2015-06-01-Euronoise/Figures/output/';
else
    outputdir = 'D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20150603-Euronoise\Figures\output-new\';
end

bFig3 = 1;
bFig4 = 0; 
bFig5 = 0; 
bFig6 = 0; 

plotOptions.I_FontSize = 24;
plotOptions.I_Width  = 10;
plotOptions.I_Height = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3
if bFig3
   
    plotOptsFig3 = plotOptions;
    plotOptsFig3.I_Width    = 16;
    plotOptsFig3.I_Height   = 18;
    plotOptsFig3.I_FontSize = 24;
    
    fnameidx = 2;

    acmode = 2;

    % F0
    h = r20150528_get_fig(acmode,1,fnameidx);
    hs = subplot(2,1,1);
    title('(a) Anechoic')
    hs = subplot(2,1,2);
    title('\Delta F0')
    xlabel('Time [s]')
    hf = Figure2paperfigureT(h,2,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    h = r20150528_get_fig(acmode,0,fnameidx);
    hs = subplot(2,1,1);
    title('(b) Reverberant')
    hs = subplot(2,1,2);
    title('\Delta F0')
    xlabel('Time [s]')
    
    hf = Figure2paperfigureT(h,2,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    % ac mode 4
    acmode = 4;

    % F0
    h = r20150528_get_fig(acmode,1,fnameidx);
    hs = subplot(2,1,1);
    title('(c) Anechoic')
    hs = subplot(2,1,2);
    xlabel('Time [s]')
    title('\Delta F0')
    hf = Figure2paperfigureT(h,2,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    h = r20150528_get_fig(acmode,0,fnameidx);
    hs = subplot(2,1,1);
    title('(d) Reverberant')
    hs = subplot(2,1,2);
    title('\Delta F0')
    xlabel('Time [s]')
    hf = Figure2paperfigureT(h,2,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4
if bFig4
    fnameidx = 3; % Loudness

    acmode = 2;

    h = r20150528_get_fig(acmode,1,fnameidx);
    title( '(a) Anechoic' );
    xlabel('Time [s]')
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    h = r20150528_get_fig(acmode,0,fnameidx);
    title( '(b) Reverberant' );
    xlabel('Time [s]')
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    % Figure 2.c, 2.d

    acmode = 4;
    h = r20150528_get_fig(acmode,1,fnameidx);
    title( '(c) Anechoic' );
    xlabel('Time [s]')
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    h = r20150528_get_fig(acmode,0,fnameidx);
    title( '(d) Reverberant' );
    xlabel('Time [s]')
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5

plotOptsFig56 = [];
plotOptsFig56 = ef(plotOptsFig56,'I_FontSize', 24);
plotOptsFig56 = ef(plotOptsFig56,'I_Width'   , 15);
title_LG = 'L_G [dB]';
acmodeYLim4 = [36    76];

if bFig5
    
    fnameidx = 4; % CB-max

    % ac-mode 4
    acmode = 4;
    
    close all
    h = r20150528_get_fig(acmode,1,fnameidx); 
    title('\Delta L_{Gmax}');
    ylim([-1.4 0.4])
    
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','Anechoic');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    close all
    h = r20150528_get_fig(acmode,0,fnameidx);
    title('\Delta L_{Gmax}');
    ylim([-1.4 0.4])
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','Reverberant');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

end

if bFig6
    fnameidx = 5; % CB-min

    % ac-mode 4
    acmode = 4;
    
    close all
    h = r20150227_get_fig(acmode,1,fnameidx); 
    title('\Delta L_{Gmin}');
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','Anechoic');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    close all
    h = r20150227_get_fig(acmode,0,fnameidx);
    title('\Delta L_{Gmin}');
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','Reverberant');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

% function opts = get_p_fig(h,name)
% 
% opts = [];
% 
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes   
%  
% if strcmp(name,'fig-fundamental-frequency-analyser-21')
% 
%     
%     xdata2 = get(dataObjs{2}, 'XData');  %data from low-level grahics objects
%     ydata2 = get(dataObjs{2}, 'YData');
%     
%     xdata1 = get(dataObjs{3}, 'XData');  %data from low-level grahics objects
%     ydata1 = get(dataObjs{3}, 'YData');
%     
%     opts.tfp  = xdata1{1};
%     opts.tfm  = xdata1{2};
%     opts.f0p = ydata1{1};
%     opts.f0m = ydata1{2};
%     
% end
% % end
% 
% function myCallbackFunction(hProp,eventData)    %#ok - hProp is unused
%    hAxes = eventData.AffectedObject;
%    tickValues = get(hAxes,'XTick');
%    newLabels = arrayfun(@(value)(sprintf('%.0f',value)), tickValues, 'UniformOutput',false);
%    set(hAxes, 'XTickLabel', newLabels);
% % end  % myCallbackFunction
% 
% function CheckPlotLims(ha)
% 
% xticks = get(gca,'XTick');
% set(gca,'XTickLabel',xticks);
% hhAxes = handle(gca);  % hAxes is the Matlab handle of our axes
% hProp = findprop(hhAxes,'XTick');  % a schema.prop object
% hListener = handle.listener(hhAxes, hProp, 'PropertyPostSet', @myCallbackFunction);
% setappdata(gca, 'XTickListener', hListener);
function r20150227_get_paper_plots
% function r20150227_get_paper_plots
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 24/02/2015
% Last update on: 24/02/2015 % Update this date manually
% Last use on   : 24/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

outputdir = '~/Documenten/Documenten-TUe/01-Text/05-Doc-TUe/lx2015-06-01-Euronoise/Figures/output/';

bFig2 = 0;
bFig3 = 0;
bFig4 = 0;
bFig5 = 1;
bFig6 = 1;

plotOptions.I_FontSize = 14;
plotOptions.I_Width  = 10;
plotOptions.I_Height = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.a, 2.b
if bFig2 == 1
    fnameidx = 1;

    acmode = 2;

    % FFT
    h = r20150227_get_fig(acmode,1,fnameidx);
    title( '(a)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    h = r20150227_get_fig(acmode,0,fnameidx);
    title( '(b)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    % Figure 2.c, 2.d

    acmode = 4;
    h = r20150227_get_fig(acmode,1,fnameidx);
    title( '(c)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    h = r20150227_get_fig(acmode,0,fnameidx);
    title( '(d)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3
if bFig3
    plotOptsFig3 = plotOptions;
    plotOptsFig3.I_Width    = 14;
    plotOptsFig3.I_FontSize = 28;
    
    fnameidx = 2;

    acmode = 2;

    % F0
    h = r20150227_get_fig(acmode,1,fnameidx);
    hs = subplot(2,1,1);
    title('(a)')
    xlabel('Time [s]')
    hs = subplot(2,1,2);
    delete(hs);
    hf = Figure2paperfigureT(h,1,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    h = r20150227_get_fig(acmode,0,fnameidx);
    hs = subplot(2,1,1);
    title('(b)')
    xlabel('Time [s]')
    hs = subplot(2,1,2);
    delete(hs);
    hf = Figure2paperfigureT(h,1,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    % ac mode 4
    acmode = 4;

    % F0
    h = r20150227_get_fig(acmode,1,fnameidx);
    hs = subplot(2,1,1);
    title('(c)')
    xlabel('Time [s]')
    hs = subplot(2,1,2);
    delete(hs);
    hf = Figure2paperfigureT(h,1,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    h = r20150227_get_fig(acmode,0,fnameidx);
    hs = subplot(2,1,1);
    title('(d)')
    xlabel('Time [s]')
    hs = subplot(2,1,2);
    delete(hs);
    hf = Figure2paperfigureT(h,1,1,plotOptsFig3);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4
if bFig4
    fnameidx = 3; % Loudness

    acmode = 2;

    % FFT
    h = r20150227_get_fig(acmode,1,fnameidx);
    title( '(a)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    h = r20150227_get_fig(acmode,0,fnameidx);
    title( '(b)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    % Figure 2.c, 2.d

    acmode = 4;
    h = r20150227_get_fig(acmode,1,fnameidx);
    title( '(c)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    h = r20150227_get_fig(acmode,0,fnameidx);
    title( '(d)' );
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5

plotOptsFig56 = [];
plotOptsFig56 = ef(plotOptsFig56,'I_FontSize', 18);
plotOptsFig56 = ef(plotOptsFig56,'I_Width'   , 15);
title_LG = 'Level L_G [dB]';
acmodeYLim2 = [11 53];
acmodeYLim4 = [31 73];

if bFig5
    fnameidx = 4; % CB-max

    acmode = 2;
    
    % option.bScale = 0;
    % option.Format = 'epsc';
    
    close all
    h = r20150227_get_fig(acmode,1,fnameidx); 
    title('(e)');
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim2);
    ht = get(ha,'Title');
    set(ht,'String','(a)');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    close all
    h = r20150227_get_fig(acmode,0,fnameidx);
    title('(f)');
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim2);
    ht = get(ha,'Title');
    set(ht,'String','(b)');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    % ac-mode 4
    acmode = 4;
    
    close all
    h = r20150227_get_fig(acmode,1,fnameidx); 
    title('(g)');
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','(c)');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    close all
    h = r20150227_get_fig(acmode,0,fnameidx);
    title('(h)');
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','(d)');
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

    acmode = 2;
    
    % option.bScale = 0;
    % option.Format = 'epsc';

    close all
    h = r20150227_get_fig(acmode,1,fnameidx); 
    title('(e)');
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim2);
    ht = get(ha,'Title');
    set(ht,'String','(a)');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    close all
    h = r20150227_get_fig(acmode,0,fnameidx);
    title('(f)');
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim2);
    ht = get(ha,'Title');
    set(ht,'String','(b)');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-rev',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');

    % ac-mode 4
    acmode = 4;
    
    close all
    h = r20150227_get_fig(acmode,1,fnameidx); 
    title('(g)');
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','(c)');
    hl = get(ha,'Xlabel');
    set(hl,'String',[]);
    hl = get(ha,'Ylabel');
    set(hl,'String',title_LG);
    
    hf = Figure2paperfigureT2(h,2,1,plotOptsFig56);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    
    close all
    h = r20150227_get_fig(acmode,0,fnameidx);
    title('(h)');
        
    figHandles = get(h,'Children');
    ha = figHandles(4);
    set(ha,'YLim',acmodeYLim4);
    ht = get(ha,'Title');
    set(ht,'String','(d)');
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
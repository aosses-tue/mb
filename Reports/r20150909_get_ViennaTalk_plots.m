function r20150909_get_ViennaTalk_plots
% function r20150909_get_ViennaTalk_plots
%
% 1. Description:
%       Produces figures as added into the Vienna Talk poster
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: r20150227_get_paper_plots
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
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
    outputdir = 'D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20150916-Vienna-talk-ORAL\Figures-new2\';
end
Mkdir(outputdir);

bFig1a = 1;
bFig1b = 0;
bFig2 = 0;

% plotOptions.I_FontSize = 14;
plotOptions.I_Width  = 10;
plotOptions.I_Height = 8;
plotOptions.I_FontSize = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1
if bFig1a
        
    plotOptions.I_Width  = 18;
    
    % ac mode 4
    fnameidx = 2;
    acmode = 4;

    % F0
    h = r20150227_get_fig(acmode,1,fnameidx);
    hs(1) = subplot(2,1,1);
    title('')
    xlabel('Time [s]')
    
    ha = gca;
    
    hs(2) = subplot(2,1,2);
    % delete(hs);
    ha(end+1) = gca;
    set(ha(2),'YLim',[-5 5]);
    
    set(ha(1),'YTick',[810:15:890])
    set(ha(2),'YTick',[-3:2:3])
    set(ha,'XTick',[2.7:0.1:3.3])
    title('')
    xlabel('Time [s]')
    
    hf = Figure2paperfigureT(h,2,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    
    legend('Recorded','Simulated');
    
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    Saveas(hf,sprintf('%s%s',outputdir,name),'fig');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1b
if bFig1b
    fnameidx = 3; % Loudness

    plotOptions.I_Legend = 'on';
    
    acmode = 4;
    h = r20150227_get_fig(acmode,1,fnameidx);
    title( '(b)' );
    xlabel('Time [s]')
    ylabel('Loudness [sone]')
    legend('recording','simulation')
    ha = gca;
    ylim([2.3 7])
    set(ha,'YTick',[2.9:0.6:7])
    set(ha,'XTick',[2.7:0.1:3.3])
    hf = Figure2paperfigureT(h,1,1,plotOptions);
    name = sprintf('%.0f-ac-mode-%.0f-ane',fnameidx,acmode);
    Saveas(hf,sprintf('%s%s',outputdir,name),'epsc');
    

end

if bFig2
    h = r20150727_get_simulated_sounds;
    name = 'simsounds';
    Saveas(h(2),sprintf('%s%s',outputdir,name),'epsc')
    
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
function [h ha] = r20150227_get_fig(acmode,bAne,fnameidx)
% function [h ha] = r20150227_get_fig(acmode,bAne,fnameidx)
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

bRev = ~bAne;
ha = [];

if ~isunix
    dir = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-06-01-Euronoise\Figures\';
else
    dir = ['~/Documenten/Documenten-TUe/01-Text/05-Doc-TUe/lx2015-06-01-Euronoise/Figures/' ];
end

if bAne
    subdir1 = sprintf('ac-%0.f-dist-ane-meas-model%s',acmode,delim);
end

if bRev
    subdir1 = sprintf('ac-%0.f-dist-rev-meas-model%s',acmode,delim);
end

switch acmode
    case 2
        fn = 424.4;
        xlimFFT     = [380 480];
        ylimFFT     = [10 55];
        ylimLoud    = [.8 3.4];
        xlimFluct   = [2.8 8.5];
        ylimFluctmax = [-1 2.5]; % diff
        ylimFluctmin = [-2 8]; % diff
        ylimF0      = [415 445];
        ylimF0diff  = [-8 5];
        
        if bRev % then reverberant 
            tlims = [3.132 4.214]; % 2 periods
            xdeltaF0 = -0.06;
        end
        if bAne
            tlims = [0.124 1.202];
            xdeltaF0 = 0.16;
        end
    case 4
        fn = 851.8;
        xlimFFT     = [750 950];
        ylimFFT     = [10 70];
        ylimLoud    = [2.5 7.5];
        xlimFluct   = [5.8 9.2];
        ylimFluctmax = [-1.5 0]; % diff
        ylimFluctmin = [3.5 5]; % diff
        ylimF0      = [800 900];
        ylimF0diff  = [-10 10];
        if bRev % then reverberant 
            tlims = [3.37 4.056]; % 2 periods
            xdeltaF0 = 0.04;
        end
        if bAne
            tlims = [2.698 3.324];
            xdeltaF0 = -0.06;
        end
end

T4 = [0.3144 0.2961*1.15];
T2 = [0.60   0.60];

fname = {   'fig-average-power-spectrum-analyser-01', ...   % 1. FFT
            'fig-fundamental-frequency-analyser-21', ...    % 2. F0
            'fig-loudness-analyser-12', ...                 % 3. Loudness
            'fig-loudness-fluctuation-max', ...             % 4.
            'fig-loudness-fluctuation-min-analyser-12' };   % 5.

        
for i = fnameidx
    figname = [dir subdir1 fname{i}];
    openfig([figname '.fig']);

    h = gcf;

    opts = get_p_fig(h,fname{i}); % we suppose data are synchronised

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(fname{i},'fig-average-power-spectrum-analyser-01')
        
        xlim(xlimFFT);
        ylim(ylimFFT);
        
        CheckPlotLims(gca);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(fname{i},'fig-fundamental-frequency-analyser-21')

        tfm = opts.tfm;
        tfp = opts.tfp;
        f0m = opts.f0m;
        f0p = opts.f0p;

        xlim(tlims);

        stPlot = [];
        
        options.trange = tlims;
        options.tanalysis = 0;
        options.delta = xdeltaF0;
        options.normfactor = fn;
        dt = tfm(2)-tfm(1);
        deltasamp = round(xdeltaF0/dt); % cant samples to delay tfp
        
        if deltasamp >= 0
            f0p = [nan(1,deltasamp) f0p(1:end-deltasamp)];
        else
            f0p = [f0p(abs(deltasamp):end) nan(1,abs(deltasamp)-1)];
        end
        
        if acmode == 2 & bRev == 1
            yoffsetF0 = -7;
            f0m = f0m+yoffsetF0;
        end
        [h ha] = Plot_fundamental_frequency(tfm,f0m,tfp,f0p,options,stPlot);
        
        set(ha(end-1),'YLim',ylimF0);
        set(ha(end)  ,'YLim',ylimF0diff);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(fname{i},'fig-loudness-analyser-12')

        % 4-dist-rev
        xlim(tlims);
        ylim(ylimLoud);

    elseif strcmp(fname{i},'fig-loudness-fluctuation-max')
        
        xlim(xlimFluct);
        ylim(ylimFluctmax);
        
    elseif strcmp(fname{i},'fig-loudness-fluctuation-min-analyser-12')
        
        xlim(xlimFluct);
        ylim(ylimFluctmin);
        
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function opts = get_p_fig(h,name)

opts = [];

axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes   
 
if strcmp(name,'fig-fundamental-frequency-analyser-21')

    
    xdata2 = get(dataObjs{2}, 'XData');  %data from low-level grahics objects
    ydata2 = get(dataObjs{2}, 'YData');
    
    xdata1 = get(dataObjs{3}, 'XData');  %data from low-level grahics objects
    ydata1 = get(dataObjs{3}, 'YData');
    
    opts.tfp  = xdata1{1};
    opts.tfm  = xdata1{2};
    opts.f0p = ydata1{1};
    opts.f0m = ydata1{2};
    
end
% end

function myCallbackFunction(hProp,eventData)    %#ok - hProp is unused
   hAxes = eventData.AffectedObject;
   tickValues = get(hAxes,'XTick');
   newLabels = arrayfun(@(value)(sprintf('%.0f',value)), tickValues, 'UniformOutput',false);
   set(hAxes, 'XTickLabel', newLabels);
% end  % myCallbackFunction

function CheckPlotLims(ha)

xticks = get(gca,'XTick');
set(gca,'XTickLabel',xticks);
hhAxes = handle(gca);  % hAxes is the Matlab handle of our axes
hProp = findprop(hhAxes,'XTick');  % a schema.prop object
hListener = handle.listener(hhAxes, hProp, 'PropertyPostSet', @myCallbackFunction);
setappdata(gca, 'XTickListener', hListener);
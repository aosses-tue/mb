function r20150227_update
% function r20150227_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 24/02/2015
% Last update on: 24/02/2015 % Update this date manually
% Last use on   : 24/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

close all
                % ac    N - anechoic    N - reverberant
                % mode  5   50  95      5   50  95
matrix4conv = [ 2       1.1  2.1  2.7     1.0  2.3     3.0; ... % measured
                2       0.9  2.0  2.6     0.9  1.9     2.9;... % modelled
                4       3.9  5.7  7.0     3.6  5.5     6.8;... % measured
                4       3.1  5.4  7.8     3.1  5.0     7.6];   % modelled
DR_ane = matrix4conv(:,4)-matrix4conv(:,2);
DR_rev = matrix4conv(:,7)-matrix4conv(:,5);

factor = 10; % 1 decimal
matrix4latex = round( factor*[matrix4conv(:,[1 2:4]) DR_ane matrix4conv(:,[5:7]) DR_rev] )/factor;

disp('%%%%%% Start copying here:')
var2latex(matrix4latex);
disp('%%%%%% End copying here: ')


acmode = 2;
bAne = 1;
bRev = ~bAne;

dir = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-06-01-Euronoise\Figures\';
if bAne
    subdir1 = sprintf('ac-%0.f-dist-ane-meas-model%s',acmode,delim);
end

if bRev
    subdir1 = sprintf('ac-%0.f-dist-rev-meas-model%s',acmode,delim);
end

switch acmode
    case 2
        fn = 424.4;
        xlimFFT     = [350 550];
        ylimFFT     = [10 60];
        ylimLoud    = [.8 3.4];
        xlimFluct   = [3 8.5];
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
        xlimFluct   = [6 9];
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

for i = 2 % 1:length(fname)
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
        % figure;
        % plot(t/T4(1),y1), hold on
        % plot(t/T4(2),y2,'r')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(fname{i},'fig-loudness-analyser-12')

        % 4-dist-rev
        xlim(tlims);
        ylim(ylimLoud);

    elseif strcmp(fname{i},'fig-loudness-fluctuation-max')
        
        xlim(xlimFluct);
        
    elseif strcmp(fname{i},'fig-loudness-fluctuation-min-analyser-12')
        
        xlim(xlimFluct);
        
    end
        
    Saveas(gcf,[figname '-f']);
    
end

if bDiary
	diary off
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
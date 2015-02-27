function [h ha] = Plot_fundamental_frequency(tfm,f0m,tfp,f0p,options,stPlot)
% function [h ha] = Plot_fundamental_frequency(tfm,f0m,tfp,f0p,options,stPlot)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 26/02/2015
% Last update on: 26/02/2015 % Update this date manually
% Last use on   : 26/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h   = [];
ha  = [];

options = Ensure_field(options,'tanalysis',0); % no-delay
% options = Ensure_field(options,'delta',0);
options = Ensure_field(options,'normfactor',1);

n = 2; 

Lf0 = min(length(tfm),length(tfp));

t1 = Do_truncate(tfm,Lf0) + options.tanalysis(1); % same time for both
t2 = t1;

options = Ensure_field(options,'trange',minmax(t1));

f0m = Do_truncate(f0m,Lf0);
f0p = Do_truncate(f0p,Lf0);

stPlot = Ensure_field(stPlot,'color',{'b-','r--'});
stPlot = Ensure_field(stPlot,'LineWidth',[1 2]);

stPlot = ef(stPlot,'Title3','(c) ');
stPlot = ef(stPlot,'Title4','(d) ');
stPlot = Ensure_field(stPlot,'bYLabel',1);

if n==2

    idx = find(t1>=options.trange(1) & t1<= options.trange(2));
    
    figure;
    subplot(n,1,1)
    plot(   t1(idx),f0m(idx), stPlot.color{1},'LineWidth',stPlot.LineWidth(1)), grid on, hold on
    plot(   t2(idx),f0p(idx), stPlot.color{2},'LineWidth',stPlot.LineWidth(2))
    grid on, hold on;
    ha(end+1) = gca;
    if stPlot.bYLabel == 1
        ylabel('Freq. [Hz]')
    end
    title(stPlot.Title3)
    ylims = get(ha(end),'YLim');
    set(ha(end),'YLim',[0.98*ylims(1) 1.02*ylims(2)]); % expand YLim in 20%            

    errorf0 = (f0p-f0m)/options.normfactor;
    if options.normfactor ~=1
        errorf0 = 100*errorf0;
    end
    subplot(n,1,2)
    plot(t1(idx), errorf0(idx), 'k')

    m = Get_mean( (abs(errorf0(idx)))' );

    legend(sprintf('avg. diff=%.2f [Hz]',m));

    grid on, hold on;
    title(stPlot.Title4)
    ha(end+1) = gca;
    xlabel('Time [s]')
    if stPlot.bYLabel == 1
        if options.normfactor == 1
            ylabel('\Delta f_0 [Hz]')
        else
            ylabel('\Delta f0/fn [%]');
        end
    end
    ylims = get(ha(end),'YLim');
    set(ha(end),'YLim',1.2*ylims); % expand YLim in 20%
    
end

linkaxes(ha,'x')
xlim(options.trange);

h(end+1) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [h ha] = Get_waveforms_and_F0_praat(f1,f2,options,stPlot)
% function [h ha] = Get_waveforms_and_F0_praat(f1,f2,options,stPlot)
%
% 1. Description:
% 
% 2. Stand-alone example:
% 
% 3. Additional info:
%   Tested cross-platform: Yes
% 
% Inspired in script: VoD_read_aligned;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file name: VoD_read_aligned
% Created on    : 21/01/2015
% Last update on: 18/02/2015 % Update this date manually
% Last use on   : 18/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h   = [];
ha  = [];

if nargin < 4
    stPlot = [];
end 

if nargin < 3
    options = [];
end 

options = Ensure_field(options,'label1','input-1');
options = Ensure_field(options,'label2','input-2');

[y1 fs] = Wavread(f1);
[y2 fs] = Wavread(f2);
L = min( length(y1), length(y2));

options = ef(options,'tanalysis',[0 (0:L-1)/fs]);

y1 = y1(1:L);
y2 = y2(1:L);
t = (0:L-1)/fs + options.tanalysis(1);

f1praat = Delete_extension(f1,'wav');
f2praat = Delete_extension(f2,'wav');

Get_F0_AC_praat([f1praat '.wav'],[f1praat '.txt.praat']);
Get_F0_AC_praat([f2praat '.wav'],[f2praat '.txt.praat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 2; 

[tfm f0m] = Get_F0_praat_from_txt( [f1praat '.txt.praat'] );
[tfp f0p] = Get_F0_praat_from_txt( [f2praat '.txt.praat'] );

Lf0 = min(length(tfm),length(tfp));

tf = Do_truncate(tfm,Lf0) + options.tanalysis(1); % same time for both

f0m = Do_truncate(f0m,Lf0);
f0p = Do_truncate(f0p,Lf0);

stPlot = Ensure_field(stPlot,'color',{'b-','r--'});
stPlot = Ensure_field(stPlot,'LineWidth',[1 2]);

stPlot = ef(stPlot,'Title1','(a) ');
stPlot = ef(stPlot,'Title2','(b) ');
if n == 4
    stPlot = ef(stPlot,'Title3','(c) ');
    stPlot = ef(stPlot,'Title4','(d) ');
else
    stPlot = ef(stPlot,'Title3','(a) ');
    stPlot = ef(stPlot,'Title4','(b) ');
end
stPlot = Ensure_field(stPlot,'bYLabel',1);

figure; 
subplot(n,1,1)
plot(t,y1,stPlot.color{1}), grid on
ha = gca;

title(stPlot.Title1)
legend(options.label1);

stPlot = rmfield(stPlot,'Title1');
if stPlot.bYLabel == 1
    ylabel('Amplitude')
end

ylims = get(ha,'YLim');
set(ha(end),'YLim',1.3*ylims); % expand YLim in 20%

subplot(n,1,2) 
plot(t,y2,'r-');, grid on
ha(end+1) = gca;
title(stPlot.Title2)
legend(options.label2);

ylims = get(ha(end),'YLim');
set(ha(end),'YLim',1.3*ylims); % expand YLim in 20%

if n == 2
    xlabel('Time [s]')
end

if stPlot.bYLabel == 1
    ylabel('Amplitude')
end

h(end+1) = gcf;

if n==2

    idx = find(tf>=options.trange(1) & tf<= options.trange(2));
    
    figure;
    subplot(n,1,1)
    plot(   tf(idx),f0m(idx), stPlot.color{1},'LineWidth',stPlot.LineWidth(1)), grid on, hold on
    plot(   tf(idx),f0p(idx), stPlot.color{2},'LineWidth',stPlot.LineWidth(2))
    grid on, hold on;
    ha(end+1) = gca;
    if stPlot.bYLabel == 1
        ylabel('Freq. [Hz]')
    end
    title(stPlot.Title3)
    ylims = get(ha(end),'YLim');
    set(ha(end),'YLim',[0.98*ylims(1) 1.02*ylims(2)]); % expand YLim in 20%            

    errorf0 = (f0p-f0m);
    subplot(n,1,2)
    plot(tf(idx), errorf0(idx), 'k')

    m = Get_mean( (abs(errorf0(idx)))' );

    legend(sprintf('avg. diff=%.2f [Hz]',m));

    grid on, hold on;
    title(stPlot.Title4)
    ha(end+1) = gca;
    xlabel('Time [s]')
    if stPlot.bYLabel == 1
        ylabel('\Delta f_0 [Hz]')
    end
    ylims = get(ha(end),'YLim');
    set(ha(end),'YLim',1.2*ylims); % expand YLim in 20%
    
end

linkaxes(ha,'x')
xlim(options.trange);


h(end+1) = gcf;
% ha(end+1) = ha(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
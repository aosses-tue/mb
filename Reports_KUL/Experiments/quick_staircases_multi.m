function outs = quick_staircases_multi(filename, opts)
% function outs = quick_staircases_multi(filename, opts)
%
% 1. Description:
%       Since 09/03/2015 this script is compatible to extract results from
%       multiprocedure experiments.
%       Valid for 1 directory at a time
% 
% 2. Stand-alone example:
%       filename = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\dau1996b-I2-3AFC-results.apr.xml';
%       opts.N4mean = 10;
%       quick_staircases_multi(filename,opts);
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2015
% Adapted from  : quick_staircases.m
% Created on    : 17/03/2015
% Last update on: 20/03/2015 % Update this date manually
% Last use on   : 20/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    opts = [];
end

if nargin < 1
    [filename p] = uigetfile({'*apr.xml','APEX result files';'*.*','All files'},'Pick up an APEX result file','D:\Documenten-TUe\02-Experiments\'); %'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\dau1996b-I2-3AFC-results.apr.xml';
    filename = [p filename];
end

opts = Ensure_field(opts,'N4mean',6);
opts = ef(opts,'bPlot',0);

N4mean  = opts.N4mean;

SRT = [];

FontSize = 14;

Stairs  = [];
ymin    = -5;
ymax    = 15;

files{1} = filename;
    
i = 1;

file = files{1};

tmp_staircase = a3adaptiveresults(file);

N = 0;
for i = 1:length(fieldnames(tmp_staircase))
    exp = sprintf('N = max(N, length(tmp_staircase.procedure%.0f));',i);
    eval(exp);
end

Stairs = [Stairs; nan(1,N)];
NumProcedures = 1; % initialisation

SRT = mean(tmp_staircase.procedure1(end-N4mean+1:end));
Stairs(end,1:length(tmp_staircase.procedure1)) = tmp_staircase.procedure1';
Stdev = std(tmp_staircase.procedure1(end-N4mean+1:end));

for k = 2:length(fieldnames(tmp_staircase)) 
    Stairs = [Stairs; nan(1,N)];
    exp = sprintf('SRT(%.0f) = mean(tmp_staircase.procedure%.0f(end-%.0f+1:end));',k,k,N4mean);
    eval(exp);
    exp = sprintf('Stairs(%.0f,1:length(tmp_staircase.procedure%.0f)) = tmp_staircase.procedure%.0f;',k,k,k);
    eval(exp)
    exp = sprintf('Stdev(%.0f) = std(tmp_staircase.procedure%.0f(end-%.0f+1:end));',k,k,N4mean);
    eval(exp);
    NumProcedures = NumProcedures + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot:
if nargout == 0 | opts.bPlot == 1
    figure
    for k = 1:NumProcedures
        subplot(NumProcedures,1,k)
        plot(1:N, Stairs(k,:),'ro','LineWidth',4), hold on

        plot([0 N],[SRT(k) SRT(k)],'k--','LineWidth',2)

        plot(1:N, Stairs(k,:)) % to plot continuous line

        hxLabel = xlabel('Trial number');
        set(hxLabel,'FontSize',FontSize)
        hyLabel = ylabel('Adapt. param (dB)');
        set(hyLabel,'FontSize',FontSize)
        hLegend = legend('staircase',['Thres = ' Num2str(SRT(k)) 'dB; std =' Num2str(Stdev(k)) ' dB']);
        set(hLegend,'Location','NorthEast')
        set(hLegend,'FontSize',FontSize)

        tmp = strsplit(filename,'.');
        tmp = strsplit(tmp{1},delim);
        name2save = name2figname( tmp{end} );

        if NumProcedures == 1
            hTitle = title(name2save);
        else
            hTitle = title( [name2save ' - procedure' num2str(k)]);
        end

        set(hTitle,'FontSize',FontSize)
        % ylim([ymin ymax])
        grid on

    end
    outs.h = gcf;
end

outs.SRT = SRT;
outs.Stdev = Stdev;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
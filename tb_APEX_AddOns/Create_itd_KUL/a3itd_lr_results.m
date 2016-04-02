function r=a3itd_lr_results(filename,logscale,type,adjustfortransducer,doplot)
% function r=a3itd_lr_results(filename,logscale,type,adjustfortransducer,doplot)
%
% 1. Description:
% 
% type: 1: least squares fit
%       2: psignifit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<2)
    logscale=0;
end

if (nargin<3)
    type=2;
end

if (nargin<4)
    adjustfortransducer=0;
end


if (nargin<5)
    doplot=1;
end

spattern='stimulus%g\t';


[labels, percent, count, totals]=a3cst2psy(filename, spattern);
%plotpsy(labels,percent, 0, [], totals);

r.data.labels=labels;
r.data.percent=percent;
r.data.count=count;
r.data.totals=totals;

%disp(sprintf('B* Total number of trials: %d', sum(totals)));
r.ntrials=sum(totals);

if (logscale)
    if (min(labels)<0)
        if (min(-labels)<0)
            warning('Can''t logscale data');
        else
            labels=-labels;
        end
    end
     set(gca, 'XScale', 'log')
end

if (adjustfortransducer)
    labels=-labels-1154;
end

if (type==1) % adaptive
    plotpsy(labels,percent);
elseif (type==2) % constant
    
    figure
    plotpsy(labels,percent);
    warning('By-passed by AO')
    % if (doplot)
    %     [cross,slope,jnd,ecross,eslope,ejnd,s]=fitpsy(labels,percent/100,totals,0,1);
    % else
    %     [cross,slope,jnd,ecross,eslope,ejnd,s]=fitpsy(labels,percent/100,totals);
    % end
    
    % r.cross=cross;
    % r.slope=slope;
    % r.jnd=jnd;
    % r.ecross=ecross;
    % r.eslope=eslope;
    % r.ejnd=ejnd;
    % r.s=s;
elseif (type==0)
    plot( [labels(1) labels(end)], [100/3 100/3], '--g');
end

xlabel('Interaural delay (us)');
ylabel('Proportion left');


if (doplot && ~iscell(filename))
    title(getshortfilename(filename));
end

% disp(sprintf('E* Total number of trials: %d  for JND=%3.1f', sum(totals), jnd));


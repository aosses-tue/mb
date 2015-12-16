function y = l0HM140_20151210_RK
% function y = l0HM140_20151210_RK
%
% 1. Description:
%       Advanced perception, assignment 4.
%       bPart1: 
%           - plots the optimal weights as a function of the variance of the visual weight
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 10/12/2015
% Last update on: 10/12/2015 
% Last use on   : 10/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bPart1 = 0;
bPart2 = 1;

if bPart1
    % Literature: Ernst and Banks 2002
    %%% Using MATLAB do the following: (1) and (2)
    
    % Source: "D:\Documenten-TUe\09-Training+activities\Advanced-perception\2015-2016-Q2\Assignment4\Processing\uitwerking assignment 4\":
    %               - uitwerking opdracht 4\script.m
    %               - uitwerking opdracht 4\Svh.m
    
    sh      = 1;        % haptic size 
    sigma2h = 1;        % haptic variance
    sv      = 2;        % visual size
    sigma2v = 0:0.1:3;  % visual variance
    [svh, sigma2vh, wv, wh] = il_Svh(sv, sigma2v, sh, sigma2h);
    
    figure;
    plot(sigma2v, [wh; wv]);
    line([sigma2h sigma2h],[0 1],'Color','r','Linestyle','--');
    text(sigma2h+0.05,0.1,'haptic variance \sigma^2_{h}');
    xlabel('visual variance \sigma^2_{v}');
    ylabel('weight values');
    legend('haptic weight','visual weight');

    figure;
    plot(sigma2v, svh);
    line([sigma2h sigma2h],[0 1.5],'Color','r','Linestyle','--');
    line([0 1],[1.5 1.5],'Color','r','Linestyle','--');
    line([0 3],[1 1],'Color','k','Linestyle',':');
    text(sigma2h+0.05,0.1,'haptic variance \sigma^2_{h}');
    xlabel('visual variance \sigma^2_{v}');
    ylabel('weighted average Svh');
    legend('haptic weight','visual weight');
    
end

if bPart2
    
    % Source: "D:\Documenten-TUe\09-Training+activities\Advanced-perception\2015-2016-Q2\Assignment4\Processing\uitwerking assignment 4\":
    %               0. - uitwerking opdracht 5\figuur1.m or figuur.m (same script)
    %               1. - uitwerking opdracht 5\script.m
    %               1. - uitwerking opdracht 5\script.m
    %               3. - uitwerking opdracht 5\posterior.m (il_posterior)
    
    %%% 0. figuur1.m
    figure;
    x=0:0.1:1;
    plot(x,[2*x; x; 0*x; -1*x; -2*x])
    hold on; plot(x,x.*x,'--')
    line([0.5 0.5],[-2 2],'Linestyle',':','Color','r')
    xlabel('x positie');
    ylabel('y positie');
    legend('i=2','i=1','i=0','i=-1','i=-2','Location','NorthWest');
    
    %%% 1. script.m
    t = 0.5;
    yhat=-2:0.01:2; % range of observed y-positions
    likelihoods=[il_lik2(yhat,-2*t,1);
                 il_lik2(yhat,-1*t,1);
                 il_lik2(yhat, 0*t,1);
                 il_lik2(yhat, 1*t,1);
                 il_lik2(yhat, 2*t,1)];
    
	figure;
    plot(yhat,likelihoods);
    xlabel('observed y-position');
    ylabel('likelihood of each goal');
    legend('i=-2','i=-1','i=0','i=1','i=2','Location','SouthWest');

    %-- Question: for what values of yhat would you infer that the goal is i = +1 
    %-- Answer: find the list indices for which the likelihood of goal i=+1 
    % is larger than the neighbouring likelihoods (i=0; i=+2). This might be
    % clearer if you look at the previous figure.    
    l1=find(il_lik2(yhat,1*t,1)>=il_lik2(yhat,2*t,1)); 
    l2=find(il_lik2(yhat,1*t,1)>=il_lik2(yhat,0*t,1)); 
    ymax=yhat(max(l1)); % intersection point (''upper'' bound)
    ymin=yhat(min(l2)); % intersection point (''lower'' bound)

    % generate output
    disp(['range of y-values that favour i=+1: (' num2str(ymin) ',' num2str(ymax) ')']);
    line([ymin ymin], [0 lik(ymin,t,1,1)],'Linestyle','--');
    line([ymax ymax], [0 lik(ymax,t,1,1)],'Linestyle','--');

    %-- Question: a posteriori
    t=0:0.001:1; % a time series

    figure;
    subplot(2,1,1);
    plot(t,[il_posterior(1,t.^2,t); il_posterior(-1,t.^2,t)]); % for a postetiori we use the measurement
    ylim([0 0.25]);
    xlabel('time t');
    ylabel('posterior probability');
    legend('i=+1','i=-1','Location','SouthWest');

    % because the priors are the same the likelihood ratio is identical to the
    % ratio of posterior probabilities
    likelihoodratio=il_lik2(t.^2,1*t,1)./il_lik2(t.^2,-1*t,1);
    l1=find(likelihoodratio>=2);
    tmin=t(min(l1));
    disp(['likelihoodratio: time after which goal i+1 is most likely: t=' num2str(tmin)]);
    posteriorratio=il_posterior(1,t.^2,t)./il_posterior(-1,t.^2,t);
    l2=find(posteriorratio>=2);
    tmin2=t(min(l2));
    disp(['posteriorratio: time after which goal i+1 is most likely: t=' num2str(tmin2)]);

    % plot it
    subplot(2,1,2);
    plot(t,[likelihoodratio' posteriorratio']);
    xlabel('time t');
    ylabel('likelihood ratio');
    line([0 1],[2 2],'LineStyle','--','Color','r');
    line([tmin tmin ],[0 likelihoodratio(min(l1)) ],'LineStyle','--','Color','r');
    line([tmin2 tmin2],[0 posteriorratio(min(l2))],'LineStyle','--','Color','r');
    text(0.1, 2.5,'threshold');
    disp('')
end

disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function [svh, sigma2vh, wv, wh]=il_Svh(sv, sigma2v, sh, sigma2h)
% function [svh, sigma2vh, wv, wh]=il_Svh(sv, sigma2v, sh, sigma2h)
%
% il_Svh returns the optimally weighted average (svh) of the visual 
% (sv, sigma2v) and haptic (sh, sigma2h) inputs.

wv      = sigma2h./(sigma2v+sigma2h);   % visual weight
wh      = 1-wv;                         % haptic weight
sigma2vh= wv.*sigma2v;                  % combined variance
svh     = wv.*sv + wh.*sh;              % weighted average

function res = il_lik(yhat,t,i,sigma)
% function res = il_lik(yhat,t,i,sigma)
%
% il_lik returns the likelihood of observing yhat at time t when the goal is
% target i and the uncertainty is sigma. The mean is y_i.

y_i = i.*t;
c   = 1/sqrt(2*pi)/sigma;
res = c*exp(-(yhat-y_i).^2/2/sigma^2);

function res = il_lik2(yhat,y_i,sigma)
% function res = il_lik2(yhat,y_i,sigma)
%
% il_lik returns the likelihood of observing yhat at time t when the goal is
% target i and the uncertainty is sigma. The mean is y_i.

c   = 1/sqrt(2*pi)/sigma;
res = c*exp(-(yhat-y_i).^2/2/sigma^2);

function res=il_posterior(i,yhat,t)
% function res=il_posterior(i,yhat,t)
%
% il_posterior returns the posterior probability of goal i for a given
% observation yhat at time t.

prior=[0 0.4 0 0.6 0]; % these are the priors for i = -2, -1, 0, 1, 2 respectively
idx = i + 3; % corresponding index. For i = -2, idx = 1; i = -1, idx = 2; etc.
sigma=1;
res=il_lik2(yhat,i*t,sigma)*prior(idx);
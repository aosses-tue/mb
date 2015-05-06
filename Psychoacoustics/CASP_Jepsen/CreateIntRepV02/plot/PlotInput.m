%--------------------------------------------------------------------------
% PlotInput.m
%--------------------------------------------------------------------------
%   plots the input signal x across time
%
% usage
%   PlotInput(time1, IntRep)
%
% input
%   time1           : time in seconds
%   IntRep          : structure containing data generated by PemoPreProc
%                     resp. CaspPreProc, browse these for details
% version 1.0
%   27/01/2013, C.T. Iben

function PlotInput(time1,IntRep)

if nargin < 2, error('not enough input arguments'); end

figH1=figure;
set(figH1,'Name','input signal','NumberTitle','off')
plot(time1,IntRep.x)
set(gca,'YLim',[1.1*min(IntRep.x) 1.1*max(IntRep.x)])
xlabel('time in s')
ylabel('amplitude')
title('input signal','Fontweight','bold')

%eof
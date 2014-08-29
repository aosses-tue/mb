function [lf le_diff outs] = ch_fluct(main_N)
% function [lf le_diff outs] = ch_fluct(main_N)
% 
% 1. Description:
%       Calculates loudness fluctuation lf for normal and hearing-impaired 
%       listeners
%
% References:
% [1] Chalupper, J. (2001) - in german - : Perzeptive Folgen von Innenohrschwerhörigkeit: 
% Modellierung, Simulation und Rehabilitation. Dissertation at the Technical
% University of Munich, Shaker Verlag.
% [2] Chalupper, J.(2000): Modellierung der Lautstärkeschwankung für Normal- und 
% Schwerhörige. Tagungsband DAGA 2000, 26. Jahrestagung Akustik, Oldenburg,
% 20.-24.3.2000, S. 254-255. 
% [3] Chalupper, J. (2007): Modeling loudness fluctuation for norm and
% hearing-impaired listeners. Proceedings of EFAS 2007, Heidelberg. in preparation 
%
% required functions: ch_korrel.m, ch_kernl2lg.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000
% Edited on     : 06/02/2007 (new version with comments and examples)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HL_ohc  = zeros(1,24);
HL_ihc  = zeros(1,24);

y_k     = ch_korrel(main_N); %cross-channel correlation
nmax    = ch_prctile(main_N,95);
nmin    = ch_prctile(main_N,5);
% figure; 
% subplot(2,1,1)
% plot([1:24], nmax, [1:24], nmin), legend('max main loudness','min main loudness')

%maximale und minimale LE's berechnen (aus 5% und 95% Lautheits-Perzentil)
le_max  = ch_kernl2lg(nmax,HL_ohc,HL_ihc); %inverse loundess transformation for normal hearing
le_min  = ch_kernl2lg(nmin,HL_ohc,HL_ihc);

% subplot(2,1,1)
% plot([1:24], le_max, [1:24], le_min), legend('max excitation pattern','min excitation pattern')

le_diff = le_max-le_min;
le_diff(find(le_diff>30))=30;  %limit to 30 dB

%account for channel correlation 
le_diff=y_k.*le_diff;

dlksum=sum(le_diff);

%transform into categorical units
% How strong (not how fast) is loudness fluctuating?
%       0=not at all, 
%       1=very weak, 
%       2=weak, 
%       3=medium, 
%       4=strong, 
%       5=very strong,
%       6=extremely strong

a=-0.19151712;
b=0.199654111;

lf=a+b.*(dlksum.^0.5);

lf(find(lf<0))=0;
lf(find(lf>6))=6;

outs.le_max = le_max;
outs.le_min = le_min;
outs.le_diff = le_diff;
outs.lf = lf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

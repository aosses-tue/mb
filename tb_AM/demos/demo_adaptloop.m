function h = demo_adaptloop
% function h = demo_adaptloop
%
% Show the effect of adaptation: This script demonstrates the effect of 
%   adaptation applied to a test signal with and without noise.
%
%   The test signal is made of a sinosoidal ramp up and down between 0 and 1.
%
%   Figure 1: Clean test signal
%
%      This figure shows the effect of adaptation on the clean test signal 
%      with and without overshoot limiting.
%
%   Figure 2: Noisy test signal
%
%      This figure shows the effect of adaptation on the noisy test signal
%      with and without overshoot limiting. Notice that in the second plot,
%      the initial spike at the beginning of the signal caused from the sharp
%      transition from complete silence to noise is magnitudes larger than
%      the values in the rest of the output.
%
%   See also: adaptloop
%
%   Url: http://amtoolbox.sourceforge.net/doc/demos/demo_adaptloop.php

% Copyright (C) 2009-2014 Peter L. Søndergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% Edited by Alejandro Osses, HTI, TU/e 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];

siglen  = 10000;
fs      = 10000;
dB_SPL  = 65; % AO
% This is the default minimum level (0 dB) of the adaptation loops. The
% loops assume that a signal is never silent, and sets all values below
% minlvl equal to minlvl. For plotting purposes, we do the same explicitly.
minlvl  = setdbspl(0);

part=siglen/10;

insig=[zeros(2*part,1);
       rampup(part);
       ones(2*part,1);
       rampdown(part);
       zeros(4*part,1)];

insig   = max(insig,minlvl);
insig   = setdbspl(insig,dB_SPL); % AO

ymindB = -80;
ymaxdB = -20;

yminMU = -250;
ymaxMU =  950;

figure;

x=(0:siglen-1)/fs;
subplot(3,1,1);
plot(x,20*log10(insig)); grid on
title('Input signal');
xlabel('time / s');
ylabel('level [dB]');
% ylabel('level / dB');
ha = gca;
ylim([ymindB ymaxdB])

subplot(3,1,2);
plot(x,adaptloop(insig,fs,0)); grid on % plot(x,adaptloop(insig,fs,0,'adt_dau'));
title('Adaptation.');
xlabel('time / s');
ylabel('level [MU]');
% ylabel('level / model units');
ha(end+1) = gca;
ylim([yminMU ymaxMU])

subplot(3,1,3);
plot(x,adaptloop(insig,fs)); grid on
title('Adaptation w. limiting.');
ylabel('level [MU]');
% ylabel('level / model units');
xlabel('time / s');
ha(end+1) = gca;
ylim([yminMU ymaxMU])

% Add a low level of noise
insig=abs(insig+0.001*randn(siglen,1));
insig=max(insig,minlvl);

h(end+1) = gcf;

figure;

subplot(3,1,1);
plot(x,20*log10(insig)); grid on
title('Input signal with added Gaussian noise.');
ylabel('level [dB]');
% ylabel('level / dB');
xlabel('time / s');
ha(end+1) = gca;
ylim([ymindB ymaxdB])

subplot(3,1,2);
plot(x,adaptloop(insig,fs,0)); grid on
title('Adaptation.');
ylabel('level [MU]');
% ylabel('level / model units');
xlabel('time / s');
ha(end+1) = gca;
ylim([yminMU ymaxMU])

subplot(3,1,3);
plot(x,adaptloop(insig,fs)); grid on
title('Adaptation w. limiting.');
ylabel('level [MU]');
% ylabel('level / model units');
xlabel('time / s');
ha(end+1) = gca;
ylim([yminMU ymaxMU])

linkaxes(ha,'x');

h(end+1) = gcf;

h(end+1) = Figure2paperfigure( [h(end-1) h(end)],3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
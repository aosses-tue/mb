function y = Learning_windows
% function y = Learning_windows
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       Learning_windows;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 12/6/2014
% Last update: 12/6/2014 % Update this date manually
% Last used: 12/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

fs = 10000;
info.fs = fs;
info.typeplot = 4;
%dur = 0.5; % seconds
K = 50;
N = 2*K;
dur = N/fs;
f0 = fs*0.10;
f1 = fs*0.16;
A0 = 1;
A1 = 0.01;
y0 = A0*Create_sin(f0,dur,fs,0); % rectangular window
y1 = A1*Create_sin(f1,dur,fs,0);
y = y0 + y1;

ymax = 50;
ymin = ymax-70;
freqfft(y,K,info);
ylim([ymin ymax])
xlims = get(gca,'XLim');
title('Rectangular')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = fs*0.105;
f1 = fs*0.16;
y0 = A0*Create_sin(f0,dur,fs,0); % rectangular window
y1 = A1*Create_sin(f1,dur,fs,0);
y = y0 + y1;
freqfft(y,K,info);
ylim([ymin      ymax])
xlim([xlims(1)  xlims(2)])
title('Rectangular')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = fs*0.1025;
f1 = fs*0.16;
y0 = A0*Create_sin(f0,dur,fs,0); % rectangular window
y1 = A1*Create_sin(f1,dur,fs,0);
y = y0 + y1;
freqfft(y,K,info);
ylim([ymin      ymax])
xlim([xlims(1)  xlims(2)])
title('Rectangular')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = fs*0.105;
f1 = fs*0.16;
y0 = A0*Create_sin(f0,dur,fs,2); % triangular window
y1 = A1*Create_sin(f1,dur,fs,2);
y = y0 + y1;
freqfft(y,K,info);
ylim([ymin      ymax])
xlim([xlims(1)  xlims(2)])
title('Triangular')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
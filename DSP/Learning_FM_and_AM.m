function y = Learning_FM_and_AM(x)
% function y = Learning_FM_and_AM(x)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/08/2014
% Last update on: 16/08/2014 % Update this date manually
% Last use on   : 16/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bDiary = 0;
Diary(mfilename,bDiary);

% I still have to calibrate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. AM

fc      = 1000; % Hz
fmod    = 70; % Hz
m       = 100;
option  = 'm'; % 
p0      = 1;

dur     = 1; % s
fs      = 44100;

[p, t]  = Create_sin_phase(fc,dur,fs);

start_phase = pi/2; % modulation starts in maximum
[pAM, pEnv] = ch_am(p,fmod,m,option,fs,start_phase);

figure;
% subplot(2,1,2)
plot( t,pEnv,t,pAM )
xlabel('time [s]')
grid on
xlim([0 2/fmod])

legend('Envelope', 'AM-signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. FM

% fc - the same
% 

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

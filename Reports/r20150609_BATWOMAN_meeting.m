function r20150609_BATWOMAN_meeting
% function r20150609_BATWOMAN_meeting
%
% 1. Description:
%       Used to create figure to be presented at Mid-term BATWOMAN meeting
%       in Graz. The last figure was manually stored with the names piano+hummer.eps
%       and piano+hummer.emf
% 
% 2. Stand-alone example:
%       r20150609_BATWOMAN_meeting;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 07/06/2015
% Last update on: 07/06/2015 % Update this date manually
% Last use on   : 07/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fpiano = 'D:\Databases\dir01-Instruments\Piano\01-Chabassier\SONS\Cd5\pressionexpe.wav';
fhummer = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\2015-02-wav-files\02-calibrated-back-as-submitted-to-Euronoise\meas-ac-2-dist-ane.wav';

[xp fsp] = Wavread(fpiano);
[xh fsh] = Wavread(fhummer);

xh = xh(1:5*fsh);

tp = (0:length(xp)-1)/fsp;
th = (0:length(xh)-1)/fsh;

xp = setdbspl(xp,60);
xh = setdbspl(xh,70);

figure;
subplot(1,2,1)
plot(tp,xp); grid on
ha = gca;
title('Piano C#5')
xlabel('Time [s]')
ylabel('Amplitude')
xlim([0 1.5])

subplot(1,2,2)
plot(th,xh); grid on
ha(end+1) = gca;
title('Hummer, acoustic mode 2')
xlabel('Time [s]')
xlim([0 0.6])

linkaxes(ha,'y');
ylim([-0.15 0.15])

h = gcf;

h(end+1) = Figure2paperfigureT2(h,1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

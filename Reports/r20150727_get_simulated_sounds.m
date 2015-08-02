function r20150727_get_simulated_sounds
% function r20150727_get_simulated_sounds
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
% Created on    : 26/07/2015
% Last update on: 26/07/2015 % Update this date manually
% Last use on   : 26/07/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fpiano = 'D:\Databases\dir01-Instruments\Piano\01-Chabassier\SONS\Cd5\pressionsimu.wav';
fhummer = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\2015-02-wav-files\02-calibrated-as-submitted-to-Euronoise\model-ac-2-dist-ane.wav';

[xp fsp] = Wavread(fpiano);
[xh fsh] = Wavread(fhummer);

xh = xh(1:5*fsh);

xp = setdbspl(xp,60);
xh = setdbspl(xh,70);

xp = [Gen_silence(50e-3,fsp); xp];
tp = (0:length(xp)-1)/fsp;
th = (0:length(xh)-1)/fsh;

tenv = tp(1:10:end);
yp = hilbert(xp(1:10:end),1000);
envp = abs(yp);
figure; plot(envp);

figure;
subplot(1,2,1)
plot(th,xh); grid on
ha = gca;
title('Hummer, acoustic mode 2')
xlabel('Time [s]')
ylabel('Pressure [Pa]')
xlim([0 0.6])

subplot(1,2,2)
plot(tp,xp); grid on, hold on
% plot(tenv,envp,'r')

ha(end+1) = gca;
title('Piano C#5')
xlabel('Time [s]')
xlim([0.05 0.34])

linkaxes(ha,'y');
ylim([-0.15 0.15])

h = gcf;

hOpts.I_KeepColor = 0;
h(end+1) = Figure2paperfigureT2(h,1,2,hOpts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ploting markers (time-onsets)

PxH     = [0.07159  0.1316];
PyH     = [0.0091   0.067];

PxP     = [0.05754  0.0636];
PyP     = [0        0.1195];

hold on;
plot(PxP,PyP,'rd','MarkerFaceColor','r','LineWidth',5)
plot([PxP(1) PxP(1)],[PyP(1) 0.15],'r--','LineWidth',2)
plot([PxP(2) PxP(2)],[PyP(2) 0.15],'r--','LineWidth',2)

disp('Select manually first subplot of the last figure (to get the axes handle into gca)')
pause()

hold on;
plot(PxH,PyH,'ro','MarkerFaceColor','r','LineWidth',5)
plot([PxH(1) PxH(1)],[PyH(1) 0.15],'r--','LineWidth',2)
plot([PxH(2) PxH(2)],[PyH(2) 0.15],'r--','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Hummer, SPLs:')
HdB = 20*log10(PyH/2e-5)
diff(HdB)

disp('Piano, SPLs:')
PdB = 20*log10(PyP/2e-5)
diff(PdB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

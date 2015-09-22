function r20150916_BATWOMAN_tutorial1
% function r20150916_BATWOMAN_tutorial1
%
% 1. Description:
%     % Out of my analysis of bPart1: 
%     %   - BP-e01 was the reference
%     %   - From higher to lower 'impact': BK-e01, BP-e01, BK-e06, BP-e06
%     %
%     % Discussed with Winfried:
%     %   - When adjusting the loudness of the audio files you are 'matching' 
%     %     their envelopes. That is visible in the plots generated in this 
%     %     script. In terms of (our) modulation filterbank terms, you matched
%     %     the lowest modulation filter (which contains envelope variations) 
%     %     forcing the rest of the filters to have the differences...
%     %   - Calibration level: approximately 78 dB SPL at the listeners' ears 
%
% 2. Stand-alone example:
%       % Run the following command in the MATLAB command window:
%       r20150916_BATWOMAN_tutorial1;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 16/09/2015
% Last update on: 16/09/2015 
% Last use on   : 22/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bPart1 = 1;
bPart2 = 0;

dir = 'D:\Downloads\audio-Winfried\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart1
    
    files = {'BP-e01.wav', 'BP-e01_loudness.wav'; ...
             'BP-e06.wav', 'BP-e06_loudness.wav'; ...
             'BK-e01.wav', 'BK-e01_loudness.wav'; ...
             'BK-e06.wav', 'BK-e06_loudness.wav'};
    
    x = [];
    xL = []; % loudness
    Chno = 1;
    % figure(i)
    for i = 1:size(files,1)
        [xtmp fs] = Wavread([dir files{i,1}]);
        x = [x xtmp(:,Chno)];
        xrms(i) = rmsdb(x(:,end))+100;
        
        [xtmp fs] = Wavread([dir files{i,2}]);
        xL = [xL xtmp(:,Chno)];
        xrmsL(i) = rmsdb(xL(:,end))+100;
        
        t = ( 1:size(x,1) )/fs; 
        figure;
        plot(t,xL(:,end)); hold on
        legend1 = sprintf('%.1f [dB]',xrmsL(end)); 
        
        plot(t,x(:,end),'r'); grid on
        legend2 = sprintf('%.1f [dB]',xrms(end));
        legend(legend1,legend2)

    end
   
end

xenv = abs(hilbert(x));
% xenv = abs(x);

%%%
f0 = 2.5; % cut-off frequency of LPF [Hz], the same as in the lowest modulation filterbank
[b a] = IRIfolp(f0,fs); % generates unity gain LPF

xenv = filter(b,a,xenv);
figure;
plot(t,xenv) % according to this plot, the amplitude envelope is quite fast
legend(files{:,1})
title(sprintf('Envelope LP with fc at %.1f [Hz]',f0))
xlabel('Time [s]'); grid on
ha = gca;

%%%
xenvL = abs(hilbert(xL));

f0 = 2.5; 
[b a] = IRIfolp(f0,fs);

xenvL = filter(b,a,xenvL);
figure;
plot(t,xenvL) % according to this plot, the amplitude envelope is quite fast (varies very rapidly)
legend(files{:,1})
title(sprintf('Loudness balanced: Envelope LP with fc at %.1f [Hz]',f0))
xlabel('Time [s]'); grid on
ha(end+1) = gca;

linkaxes(ha,'xy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart2

fs = 2000;
tf = 3; % seconds
t = 0:1/fs:tf; % time vector
f0 = 100; % Hz
f1 = 200; % Hz
method = 'quadratic';
t1 = 1; % seconds
x = chirp(t,f0,t1,f1,method); % % Start @ f1, cross f2 Hz at t1 

figure;
spectrogram(x,256,240,256,fs);

end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

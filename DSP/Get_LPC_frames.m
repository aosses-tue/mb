function Get_LPC_frames(x,fs,N,numframes)
% function Get_LPC_frames(x,fs,N,numframes)
%
% 1. Description:
%       Plot formants using LPC. For the voice of the dragon, ideally put as
%       input 1 or 2 cycles.
% 
% 2. Stand-alone example:
%       filename= 'D:\Output\tmp-Audio\meas-ac-mode-5.wav';
%       [x fs]  = Wavread(filename);
%       N       = 4096;
%       nframes = 10;
%       Get_LPC_frames(x,fs,N,nframes);
%
%       filename= '~/Documenten/Databases/dir01-Instruments/Voice-of-dragon/03-Wav-files-predicted/04-Wav-files-calibrated-44.1kHz/modus-1-v_2filt.wav';
%       [x fs]  = Wavread(filename);
%       N       = 4096;
%       nframes = 10;
%       Get_LPC_frames(x,fs,N,nframes);
%
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 09/09/2014
% Last update on: 09/09/2014 % Update this date manually
% Last use on   : 22/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    N = 4096;
end

if nargin < 4
    numframes = 10;
end

if fs ~= 10000
    x=resample(x,10000,fs);
    fs=10000;
end

xframed = buffer(x,floor(length(x)/numframes),0);
dt = size(xframed,1)/fs;

h = [];
ha = [];

for i = 1:numframes
    
    [xx out] = Get_LPC(xframed(:,i),fs,N);
    figure;
    plot(out.f,out.h_dB - max(out.h_dB)), hold on;
    h(end+1)    = gcf;
    ha(end+1)   = gca;
    if i ~= 1 & i ~= numframes
        title([sprintf('t%.0f = ti + %.0f',i,i-1) '\Delta t'])
    elseif i == 1
        title(['ti = t0; \Delta t = ' sprintf('%.3f [s]',dt)])
    elseif i == numframes
        title(['tf = ti + ' sprintf('%.3f [s]',dt*(i-1))]);
    end
    
    if i == numframes/2 | i == numframes
        xlabel('Frequency [Hz]')
    end
    
    if i == floor( numframes/4 )+1
        ylabel('Gain [dB], normalised to 0 at f1')
    end
end

linkaxes(ha,'xy');

hM = ImageSetup; 
hM.I_Matrix      = [numframes/2,numframes/2];
hM.I_FontSize    = 10; 
hM.I_FontName    = 'Arial'; 
hM.I_Width       = 8;
hM.I_Height      = 8;
hM.I_TitleInAxis = 1;
hM.I_Space       = [0.01,0.01];

hM.I_Ylim   = [-50,0]; % Uncomment for fixing the limits in the y-axis
hM.I_Xlim    = [min(out.f) max(out.f)];

% h.I_Xlim = [0,5];
hM.I_Grid = 'on'; 
hM.I_KeepColor = 0; 
hM.I_Handles = h;
hM.prepareAllFigures;
hM.arrayAddedHandles = 1;
add2ArraySubplotVer(hM);

close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

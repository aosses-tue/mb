function demo_dau1997b_one_sound(insig,fs,options)
% function demo_dau1997b_one_sound(insig,fs,options)
%
% 1. Description:
%       Recreates simulations as presented in Dau1996b. If stimuli are not
%       found, then they are along this script generated.
% 
% 2. Stand-alone example:
%       dir = 'D:\Databases\dir01-Instruments\Piano\04-Antoine\F1\';
%       f1 = [dir 'NS19-F1.wav']; 
%       [insig fs] = Wavread(f1);
%       options.fc2plot = 90;
%       insig = insig(1:3*fs); % 3 seconds
%       demo_dau1997b_one_sound(insig,fs,options);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 01/04/2015
% Last update on: 01/04/2015 % Update this date manually
% Last use on   : 01/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    close all
    options = [];
end

options = Ensure_field(options, 'fc2plot',5000);
options = Ensure_field(options, 'bSave', 0);
options = Ensure_field(options, 'Title','');
options = Ensure_field(options, 'output_dir', Get_TUe_paths('outputs'));
paths.outputs   = options.output_dir;
fc2plot = options.fc2plot;

h = []; % we initialise handle for Figures

bListen = 0;
if nargin < 2
    fs = 44100;
end

if nargin == 0 
    %% Dau1997b, Fig. 7:
    % Common parameters:
    BW = 3; %[3 31 314];
    mdept = 100;%d2m(d,'dau');

    fc_tone = 5000;
    SPL = 65;
    dur = 3;
    fm  = 20;
    insig_NBN   = AM_random_noise_BW(fc_tone,BW,SPL,dur,fs,fm,mdept);
    insig   = [Gen_silence(0.5,fs); insig_NBN; Gen_silence(0.5,fs)];
    fc2plot = fc_tone;
    
end

t = ( 0:length(insig)-1 )/fs;

if bListen == 1
    sound(insig,fs);
end

%%

[outsig1 , fc ,fcm] = dau1997preproc(insig,fs);
fc2look = abs(fc-fc2plot);
fcmin   = min(fc2look);
idx     = find(fc2look == fcmin);
outsig1 = outsig1{idx};
fc2plot = fc(idx);

NCh_mod = size(outsig1,2);
fmod_num= 1:NCh_mod;
fmod    = fcm(fmod_num);

Zmin = min(min(outsig1(:,fmod_num)));
Zmax = max(max(outsig1(:,fmod_num)));

opts.step1  = 1;
opts.step2  = 1;
opts.Title  = sprintf('%s, fc=%.0f [Hz]',options.Title,fc2plot); 
opts.bPlot3D= 0;
opts.bPlot2D= ~opts.bPlot3D;
opts.XLabel = 'Time [s]';
opts.YLabel = 'Modulation frequency [Hz]';
opts.ZLabel = 'Amplitude [MU]';

if opts.bPlot2D
    if NCh_mod <= 4
        opts.I_Matrix = [1,4];
    elseif NCh_mod <= 8
        opts.I_Matrix = [2,4];
    else
        opts.I_Matrix = [3,4];
    end
    opts.I_Width  = opts.I_Matrix(2)*4;
    opts.I_Height = opts.I_Matrix(1)*4;
end
                                                                                                                                                                                                       
figure;
Mesh(t,fmod,transpose(outsig1(:,fmod_num)),opts);
% surfc(repmat(t,length(fmod_num),1),repmat(fmod_num',1,length(t)),transpose(outsig1(:,fmod_num)));
if opts.bPlot3D
    xlabel('Time [s]')
    ylabel('Modulation filter');
    zlabel('[MU]')
    ha = gca;
    set(ha,'ZLim',[Zmin Zmax])
    set(ha,'View',[50 35])
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

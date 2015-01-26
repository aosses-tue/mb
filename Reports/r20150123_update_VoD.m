function y = r20150123_update_VoD
% function y = r20150123_update_VoD
%
% 1. Description:
%       Data presented in report generated in GUI PsySoundControl.fig. Use
%       the screenshots saved inside the folder of each condition to see
%       the details of each configuration.
% 
%       In parallel, I started to generate the test battery for Fluctuation
%       Strength (nothing written about this by my side)
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/01/2015
% Last update on: 23/01/2015 % Update this date manually
% Last use on   : 23/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all

bPart1 = 0;
bPartFluct = 0;
bPartRough = 1;

if bPart1
% directory = '/home/alejandro/Documenten/Databases/dir01-Instruments/Mech-Hummer/';
directory = '~/Documenten/Databases/dir01-Instruments/new-ac-mode2-close+far/';

f1 = [directory 'model-ac-mode-2.wav'];
f2 = [directory 'meas-ac-mode-2.wav']; % file already aligned
f3 = [directory 'model-ac-mode-2-far.wav'];
f4 = [directory 'meas-ac-mode-2-far.wav']; 

acmode  = 5;

% model
[x fs]  = Wavread(f1);
x2      = Wavread(f2);
x3      = Wavread(f3);
x4      = Wavread(f4);

t       = ( 1:length(x) )/fs;
misc    = Get_VoD_params(0);
T       = misc.Tmodel(acmode-1); % misc.ti_model(acmode-1);

% 8 ac mode 2
% 3 ac mode 5
ti = 3*T; % starting from 8th period
tf = ti + 4*T; % ending 4 periods later

tis = ti*fs
tfs = tf*fs

sprintf('ti = %.0f; tf = %.0f',tis,tfs) % Put this in PsySoundControl.fig

% measured
offset = 0.04; % just to plot

figure;
plot(t,x +3*offset), hold on, grid on
plot(t,x2+1*offset,'r')

plot(t,x3-1*offset)
plot(t,x4-3*offset,'r')

xlim([ti tf])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cal file: white noise:

dir_where_ref = Get_TUe_paths('db_fastl2007');
dir_output  = '~/Documenten/MATLAB/outputs/tones/';

file    = [dir_where_ref 'track_03.wav'];

if bPartFluct
    
    % For Figure 10.a (Fastl 2007)
    t       = 8;
    
    [insig fs] = Wavread(file);
    options.callevel = 60;
    
    fmod    = [1 2 4 8 16 32];
    m       = 100;
    option  = 'm';
    % lvlref  = 60;
    % lvl     = 60;
    
    start_phase = -pi/2; % pi/2;
    
    for i = 1:length(fmod)
        [y env] = ch_am(insig,fmod(i),m,option,fs,start_phase);
        Wavwrite( y(1:t*fs),fs, sprintf('%sbbn_AM_m_%s_fmod_%sHz',dir_output,Num2str(m,3),Num2str(fmod(i),3)) );
    end
end

if bPartRough
    
    % For Figure 11.4.a (Fastl2007)
    t       = 1;
    
    [insig fs] = Wavread(file);
    % options.callevel = 60;
    
    fmod    = [30 50 70 90 110];
    m       = 100;
    option  = 'm';
    lvlref  = 60;
    lvl     = 60;
    
    start_phase = -pi/2; % pi/2;
    
    for i = 1:length(fmod)
        [y env] = ch_am(insig,fmod(i),m,option,fs,start_phase);
        Wavwrite( y(1:t*fs),fs, sprintf('%sbbn_AM_m_%s_fmod_%sHz',dir_output,Num2str(m,3),Num2str(fmod(i),3)) );
    end
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
